"""
Infer ancestral sequences based on a tree.

The ancestral sequences are inferred using `TreeTime <https://academic.oup.com/ve/article/4/1/vex042/4794731>`_.
Each internal node gets assigned a nucleotide sequence that maximizes a
likelihood on the tree given its descendants and its parent node.
Each node then gets assigned a list of nucleotide mutations for any position
that has a mismatch between its own sequence and its parent's sequence.
The node sequences and mutations are output to a node-data JSON file.

If amino acid options are provided, the ancestral amino acid sequences for each
requested gene are inferred with the same method as the nucleotide sequences described above.
The inferred amino acid mutations will be included in the output node-data JSON
file, with the format equivalent to the output of `augur translate`.

The nucleotide and amino acid sequences are inferred separately in this command,
which can potentially result in mismatches between the nucleotide and amino
acid mutations. If you want amino acid mutations based on the inferred
nucleotide sequences, please use `augur translate`.

.. note::

    The mutation positions in the node-data JSON are one-based.
"""
from augur.errors import AugurError
import sys
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .utils import read_tree, InvalidTreeError, write_json, get_json_name
from treetime.vcf_utils import read_vcf, write_vcf
from collections import defaultdict

def ancestral_sequence_inference(tree=None, aln=None, ref=None, infer_gtr=True,
                                 marginal=False, fill_overhangs=True, infer_tips=False,
                                 alphabet='nuc'):
    """infer ancestral sequences using TreeTime

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree or str
        tree or filename of tree
    aln : Bio.Align.MultipleSeqAlignment or str
        alignment or filename of alignment
    ref : str, optional
        reference sequence to pass to TreeTime's TreeAnc class
    infer_gtr : bool, optional
        Description
    marginal : bool, optional
        Description
    fill_overhangs : bool
       In some cases, the missing data on both ends of the alignment is
       filled with the gap character ('-'). If set to True, these end-gaps are
       converted to "ambiguous" characters ('N' for nucleotides, 'X' for
       aminoacids). Otherwise, the alignment is treated as-is
    infer_tips : bool
        Since v0.7, TreeTime does not reconstruct tip states by default.
        This is only relevant when tip-state are not exactly specified, e.g. via
        characters that signify ambiguous states. To replace those with the
        most-likely state, set infer_tips=True
    alphabet : str
        alphabet to use for ancestral sequence inference. Default is the nucleotide
        alphabet that included a gap character 'nuc'. Alternative is `aa` for amino
        acids.

    Returns
    -------
    treetime.TreeAnc
        treetime.TreeAnc instance
    """

    from treetime import TreeAnc

    tt = TreeAnc(tree=tree, aln=aln, ref=ref, gtr='JC69', alphabet=alphabet,
                 fill_overhangs=fill_overhangs, verbose=1)

    # convert marginal (from args.inference) from 'joint' or 'marginal' to True or False
    bool_marginal = (marginal == "marginal")

    # only infer ancestral sequences, leave branch length untouched
    tt.infer_ancestral_sequences(infer_gtr=infer_gtr, marginal=bool_marginal,
                                 reconstruct_tip_states=infer_tips)

    return tt

def collect_mutations_and_sequences(tt, infer_tips=False, full_sequences=False, character_map=None, is_vcf=False):
    """iterates of the tree and produces dictionaries with
    mutations and sequences for each node.

    Parameters
    ----------
    tt : treetime.TreeTime
        instance of treetime with valid ancestral reconstruction
    infer_tips : bool, optional
        if true, request the reconstructed tip sequences from treetime, otherwise retain input ambiguities
    full_sequences : bool, optional
        if true, add the full sequences
    character_map : None, optional
        optional dictionary to map characters to a custom set.

    Returns
    -------
    dict
        dictionary of mutations and sequences
    """
    if character_map is None:
        cm = lambda x:x
    else:
        cm = lambda x: character_map.get(x, x)

    data = defaultdict(dict)
    inc = 1 # convert python numbering to start-at-1
    for n in tt.tree.find_clades():
        data[n.name]['muts'] = [a+str(int(pos)+inc)+cm(d)
                                for a,pos,d in n.mutations]

    if is_vcf:
        mask = np.zeros(tt.sequence_length, dtype=bool)
    else:
        # Identify sites for which every terminal sequence is ambiguous.
        # These sites will be masked to prevent rounding errors in the
        # maximum likelihood inference from assigning an arbitrary
        # nucleotide to sites at internal nodes.
        ambiguous_count = np.zeros(tt.sequence_length, dtype=int)
        for n in tt.tree.get_terminals():
            ambiguous_count += np.array(tt.sequence(n,reconstructed=False, as_string=False)==tt.gtr.ambiguous, dtype=int)
        mask = ambiguous_count==len(tt.tree.get_terminals())

    if full_sequences:
        for n in tt.tree.find_clades():
            try:
                tmp = tt.sequence(n,reconstructed=infer_tips, as_string=False)
                tmp[mask] = tt.gtr.ambiguous
                data[n.name]['sequence'] = "".join(tmp)
            except:
                print("No sequence available for node ",n.name)

    return {"nodes": data, "mask": mask}

def run_ancestral(T, aln, root_sequence=None, is_vcf=False, full_sequences=False, fill_overhangs=False,
                  infer_ambiguous=False, marginal=False, alphabet='nuc'):
    tt = ancestral_sequence_inference(tree=T, aln=aln, ref=root_sequence if is_vcf else None, marginal=marginal,
                                      fill_overhangs = fill_overhangs, alphabet=alphabet,
                                      infer_tips = infer_ambiguous)

    character_map = {}
    for x in tt.gtr.profile_map:
        if tt.gtr.profile_map[x].sum()==tt.gtr.n_states:
            # TreeTime treats all characters that are not valid IUPAC nucleotide chars as fully ambiguous
            # To clean up auspice output, we map all those to 'N'
            character_map[x] = 'N'
        else:
            character_map[x] = x
    # add reference sequence to json structure. This is the sequence with
    # respect to which mutations on the tree are defined.
    if root_sequence:
        root_seq = root_sequence
    else:
        root_seq = tt.sequence(T.root, as_string=True)

    mutations = collect_mutations_and_sequences(tt, full_sequences=full_sequences,
                          infer_tips=infer_ambiguous, character_map=character_map, is_vcf=is_vcf)
    if root_sequence:
        for pos, (root_state, tree_state) in enumerate(zip(root_sequence, tt.sequence(tt.tree.root, reconstructed=infer_ambiguous, as_string=True))):
            if root_state != tree_state:
                mutations['nodes'][tt.tree.root.name]['muts'].append(f"{root_state}{pos+1}{tree_state}")

    return {'tt': tt,
            'root_seq': root_seq,
            'mutations': mutations}


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("ancestral", help=__doc__)

    input_group = parser.add_argument_group(
        "inputs",
        "Tree and sequences to use for ancestral reconstruction"
    )
    input_group.add_argument('--tree', '-t', required=True, help="prebuilt Newick")
    input_group.add_argument('--alignment', '-a', help="alignment in FASTA or VCF format")
    input_group.add_argument('--vcf-reference', type=str, help='FASTA file of the sequence the VCF was mapped to (only used if a VCF is provided as the alignment)')
    input_group.add_argument('--root-sequence', type=str, help='FASTA/genbank file of the sequence that is used as root for mutation calling.'
                        ' Differences between this sequence and the inferred root will be reported as mutations on the root branch.')

    global_options_group = parser.add_argument_group(
        "global options",
        "Options to configure reconstruction of both nucleotide and amino acid sequences"
    )
    global_options_group.add_argument('--inference', default='joint', choices=["joint", "marginal"],
                                    help="calculate joint or marginal maximum likelihood ancestral sequence states")

    nucleotide_options_group = parser.add_argument_group(
        "nucleotide options",
        "Options to configure reconstruction of ancestral nucleotide sequences"
    )
    ambiguous = nucleotide_options_group.add_mutually_exclusive_group()
    ambiguous.add_argument('--keep-ambiguous', action="store_true",
                                help='do not infer nucleotides at ambiguous (N) sites on tip sequences (leave as N).')
    ambiguous.add_argument('--infer-ambiguous', action="store_true", default=True,
                                help='infer nucleotides at ambiguous (N,W,R,..) sites on tip sequences and replace with most likely state.')
    nucleotide_options_group.add_argument('--keep-overhangs', action="store_true", default=False,
                                help='do not infer nucleotides for gaps (-) on either side of the alignment')

    amino_acid_options_group = parser.add_argument_group(
        "amino acid options",
        "Options to configure reconstruction of ancestral amino acid sequences. All arguments are required for ancestral amino acid sequence reconstruction."
    )
    amino_acid_options_group.add_argument('--annotation',
                        help='GenBank or GFF file containing the annotation')
    amino_acid_options_group.add_argument('--genes', nargs='+', help="genes to translate (list or file containing list)")
    amino_acid_options_group.add_argument('--translations', type=str, help="translated alignments for each CDS/Gene. "
                           "Currently only supported for FASTA-input. Specify the file name via a "
                           "template like 'aa_sequences_%%GENE.fasta' where %%GENE will be replaced "
                           "by the gene name.")

    output_group = parser.add_argument_group(
        "outputs",
        "Outputs supported for reconstructed ancestral sequences"
    )
    output_group.add_argument('--output-node-data', type=str, help='name of JSON file to save mutations and ancestral sequences to')
    output_group.add_argument('--output-sequences', type=str, help='name of FASTA file to save ancestral nucleotide sequences to (FASTA alignments only)')
    output_group.add_argument('--output-translations', type=str, help="name of the FASTA file(s) to save ancestral amino acid sequences to. "
                        "Specify the file name via a template like 'ancestral_aa_sequences_%%GENE.fasta' where %%GENE will be replaced by"
                        "the gene name.")
    output_group.add_argument('--output-vcf', type=str, help='name of output VCF file which will include ancestral seqs')

    return parser

def run(args):
    # Validate arguments.
    aa_arguments = (args.annotation, args.genes, args.translations)
    if any(aa_arguments) and not all(aa_arguments):
        raise AugurError("For amino acid sequence reconstruction, you must provide an annotation file, a list of genes, and a template path to amino acid sequences.")

    # check alignment type, set flags, read in if VCF
    is_vcf = any([args.alignment.lower().endswith(x) for x in ['.vcf', '.vcf.gz']])
    ref = None

    try:
        T = read_tree(args.tree)
    except (FileNotFoundError, InvalidTreeError) as error:
        print("ERROR: %s" % error, file=sys.stderr)
        return 1

    import numpy as np
    missing_internal_node_names = [n.name is None for n in T.get_nonterminals()]
    if np.all(missing_internal_node_names):
        print("\n*** WARNING: Tree has no internal node names!", file=sys.stderr)
        print("*** Without internal node names, ancestral sequences can't be linked up to the correct node later.", file=sys.stderr)
        print("*** If you want to use 'augur export' or `augur translate` later, re-run this command with the output of 'augur refine'.", file=sys.stderr)
        print("*** If you haven't run 'augur refine', you can add node names to your tree by running:", file=sys.stderr)
        print("*** augur refine --tree %s --output-tree <filename>.nwk"%(args.tree) , file=sys.stderr)
        print("*** And use <filename>.nwk as the tree when running 'ancestral', 'translate', and 'traits'", file=sys.stderr)

    if is_vcf:
        if not args.vcf_reference:
            print("ERROR: a reference Fasta is required with VCF-format alignments", file=sys.stderr)
            return 1

        compress_seq = read_vcf(args.alignment, args.vcf_reference)
        aln = compress_seq['sequences']
        ref = compress_seq['reference']
    else:
        aln = args.alignment
        ref = None
        if args.root_sequence:
            for fmt in ['fasta', 'genbank']:
                try:
                    ref = str(SeqIO.read(args.root_sequence, fmt).seq)
                    break
                except:
                    pass
            if ref is None:
                print(f"ERROR: could not read root sequence from {args.root_sequence}", file=sys.stderr)
                return 1

    import treetime
    print("\nInferred ancestral sequence states using TreeTime:"
          "\n\tSagulenko et al. TreeTime: Maximum-likelihood phylodynamic analysis"
          "\n\tVirus Evolution, vol 4, https://academic.oup.com/ve/article/4/1/vex042/4794731\n")

    print(f"augur ancestral is using TreeTime version {treetime.version}")

    # Infer ambiguous bases if the user has requested that we infer them (either
    # explicitly or by default) and the user has not explicitly requested that
    # we keep them.
    infer_ambiguous = args.infer_ambiguous and not args.keep_ambiguous
    full_sequences = not is_vcf
    nuc_result = run_ancestral(T, aln, root_sequence=ref if ref else None, is_vcf=is_vcf, fill_overhangs=not args.keep_overhangs,
                               full_sequences=full_sequences, marginal=args.inference, infer_ambiguous=infer_ambiguous, alphabet='nuc')
    anc_seqs = nuc_result['mutations']
    anc_seqs['reference'] = {'nuc': nuc_result['root_seq']}

    if anc_seqs.get("mask") is not None:
        anc_seqs["mask"] = "".join(['1' if x else '0' for x in anc_seqs["mask"]])

    anc_seqs['annotations'] = {'nuc': {'start': 1, 'end': len(anc_seqs['reference']['nuc']),
                                       'strand': '+', 'type': 'source'}}

    if not is_vcf and args.genes:
        from .utils import load_features
        ## load features; only requested features if genes given
        features = load_features(args.annotation, args.genes)
        if features is None:
            print("ERROR: could not read features of reference sequence file")
            return 1
        print("Read in {} features from reference sequence file".format(len(features)))
        for gene in args.genes:
            print(f"Processing gene: {gene}")
            fname = args.translations.replace("%GENE", gene)
            feat = features[gene]
            root_seq = str(feat.extract(Seq(ref)).translate()) if ref else None

            aa_result = run_ancestral(T, fname, root_sequence=root_seq, is_vcf=is_vcf, fill_overhangs=not args.keep_overhangs,
                                        marginal=args.inference, infer_ambiguous=infer_ambiguous, alphabet='aa')
            if aa_result['tt'].data.full_length*3 != len(feat):
                print(f"ERROR: length of translated alignment for {gene} does not match length of reference feature."
                       " Please make sure that the annotation matches the translations.")
                return 1

            for key, node in anc_seqs['nodes'].items():
                if 'aa_muts' not in node: node['aa_muts'] = {}
                node['aa_muts'][gene] = aa_result['mutations']['nodes'][key]['muts']

                # Add amino acid sequences to the root node of the tree.
                if key == T.root.name:
                    if "aa_sequences" not in node:
                        node["aa_sequences"] = {}

                    node["aa_sequences"][gene] = aa_result['tt'].sequence(T.root, as_string=True, reconstructed=True)

            anc_seqs['reference'][gene] = aa_result['root_seq']
            # FIXME: Note that this is calculating the end of the CDS as 3*length of translation
            # this is equivalent to the annotation for single segment CDS, but not for cds
            # with splicing and slippage. But auspice can't handle the latter at the moment.
            anc_seqs['annotations'][gene] = {'seqid':args.annotation,
                                             'type':feat.type,
                                             'start': int(feat.location.start)+1,
                                             'end': int(feat.location.start) + 3*len(anc_seqs['reference'][gene]),
                                             'strand': {+1:'+', -1:'-', 0:'?', None:None}[feat.location.strand]}

            # Save ancestral amino acid sequences to FASTA.
            if args.output_translations:
                with open(args.output_translations.replace("%GENE", gene), "w", encoding="utf-8") as oh:
                    for node in aa_result["tt"].tree.find_clades():
                        oh.write(f">{node.name}\n{aa_result['tt'].sequence(node, as_string=True, reconstructed=True)}\n")

    out_name = get_json_name(args, '.'.join(args.alignment.split('.')[:-1]) + '_mutations.json')
    write_json(anc_seqs, out_name)
    print("ancestral mutations written to", out_name, file=sys.stdout)

    if args.output_sequences:
        if args.output_vcf:
            # TODO: This should be an error and we should check for this
            # unsupported combination of arguments at the beginning of the
            # script to avoid wasting time for users.
            print("WARNING: augur only supports sequence output for FASTA alignments and not for VCFs.", file=sys.stderr)
        else:
            records = [
                SeqRecord(Seq(node_data["sequence"]), id=node_name, description="")
                for node_name, node_data in anc_seqs["nodes"].items()
            ]
            SeqIO.write(records, args.output_sequences, "fasta")
            print("ancestral sequences FASTA written to", args.output_sequences, file=sys.stdout)

    # If VCF, output VCF including new ancestral seqs
    if is_vcf:
        if args.output_vcf:
            vcf_fname = args.output_vcf
        else:
            vcf_fname = '.'.join(args.alignment.split('.')[:-1]) + '.vcf'
        write_vcf(nuc_result['tt'].get_tree_dict(keep_var_ambigs=True), vcf_fname)
        print("ancestral sequences as vcf-file written to",vcf_fname, file=sys.stdout)

    return 0
