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
from .utils import parse_genes_argument, read_tree, InvalidTreeError, write_json, get_json_name, \
    genome_features_to_auspice_annotation
from .io.file import open_file
from .io.vcf import is_vcf as is_filename_vcf
from treetime.vcf_utils import read_vcf, write_vcf
from collections import defaultdict
from .types import ValidationMode
from .util_support.node_data_file import NodeDataObject
from .export_v2 import validation_mode_help_message

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

def create_mask(is_vcf, tt, reference_sequence, aln):
    """
    Identify sites for which every terminal sequence is ambiguous. These sites
    will be masked to prevent rounding errors in the maximum likelihood
    inference from assigning an arbitrary nucleotide to sites at internal nodes.

    Parameters
    ----------
    is_vcf : bool
    tt : treetime.TreeTime
        instance of treetime with valid ancestral reconstruction. Unused if is_vcf.
    reference_sequence : str
        only used if is_vcf
    aln : dict
        describes variation (relative to reference) per sample. Only used if is_vcf.
        
    Returns
    -------
    numpy.ndarray(bool)
    """
    num_tips = len(tt.tree.get_terminals())
    ambiguous_count = np.zeros(tt.sequence_length, dtype=int)
    if is_vcf:
        variable_sites = set()
        # VCF ambiguous positions come in two forms:
        # Firstly, if VCF defines a "N" ALT and assigns it to every sample:
        for sample_data in aln.values():
            ambig_positions = []
            for pos, alt in sample_data.items():
                variable_sites.add(pos)
                if alt=='N':
                    ambig_positions.append(pos)
            # ambig_positions = [pos for pos, alt in sample_data.items() if alt=='N']
            np.add.at(ambiguous_count, ambig_positions, 1)       
        # Secondly, if the VCF defines no mutations but the ref is "N":
        no_info_sites = np.array(list(reference_sequence)) == 'N'
        no_info_sites[list(variable_sites)] = False
        # and the mask is the union of these two forms
        mask = ambiguous_count==num_tips
        mask[no_info_sites] = True
    else:
        for n in tt.tree.get_terminals():
            ambiguous_count += np.array(tt.sequence(n,reconstructed=False, as_string=False)==tt.gtr.ambiguous, dtype=int)
        mask = ambiguous_count==num_tips
    return mask

def collect_mutations(tt, mask, character_map=None, reference_sequence=None, infer_ambiguous=False):
    """iterates of the tree and produces dictionaries with
    mutations and sequences for each node.

    If a reference sequence is provided then mutations can be collected for the
    root node. Masked positions at the root-node will be treated specially: if
    we infer ambiguity, then we report no mutations (i.e. we assume the
    reference base holds), otherwise we'll report a mutation from the <ref> to
    "N".

    Parameters
    ----------
    tt : treetime.TreeTime
        instance of treetime with valid ancestral reconstruction
    mask : numpy.ndarray(bool)
    character_map : None, optional
        optional dictionary to map characters to a custom set.
    reference_sequence : str, optional

    Returns
    -------
    dict
        dict -> <node_name> -> [mut, mut, ...] where mut is a string in the form
        <from><1-based-pos><to>
    """

    if character_map is None:
        cm = lambda x:x
    else:
        cm = lambda x: character_map.get(x, x)

    data = {}
    inc = 1 # convert python numbering to start-at-1

    # Note that for mutations reported across the tree we don't have to consider
    # the mask, because while sites which are all Ns may have an inferred base,
    # there will be no variablity and thus no mutations.  
    for n in tt.tree.find_clades():
        data[n.name] = [a+str(int(pos)+inc)+cm(d)
                        for a,pos,d in n.mutations]

    if reference_sequence:
        data[tt.tree.root.name] = []
        for pos, (root_state, tree_state) in enumerate(zip(reference_sequence, tt.sequence(tt.tree.root, reconstructed=infer_ambiguous, as_string=True))):
            if mask[pos] and infer_ambiguous:
                continue
            if root_state != tree_state:
                data[tt.tree.root.name].append(f"{root_state}{pos+1}{tree_state}")

    return data
        
def collect_sequences(tt, mask, reference_sequence=None, infer_ambiguous=False):
    """
    Create a full sequence for every node on the tree. Masked positions will
    have the reference base if we are inferring ambiguity, or the ambiguous
    character 'N'.

    Parameters
    ----------
    tt : treetime.TreeTime
        instance of treetime with valid ancestral reconstruction
    mask : numpy.ndarray(bool)
        Mask these positions by changing them to the ambiguous nucleotide
    reference_sequence : str or None 
    infer_ambiguous : bool, optional
        if true, request the reconstructed sequences from treetime, otherwise retain input ambiguities

    Returns
    -------
    dict
        dict -> <node_name> -> sequence_string
    """
    sequences = {}

    ref_mask = None

    if reference_sequence and infer_ambiguous:
        ref_mask = np.array(list(reference_sequence))[mask]

    for n in tt.tree.find_clades():
        try:
            tmp = tt.sequence(n,reconstructed=infer_ambiguous, as_string=False)
            if ref_mask is not None:
                tmp[mask] = ref_mask
            else:
                tmp[mask] = tt.gtr.ambiguous
            sequences[n.name] = "".join(tmp)
        except:
            print("No sequence available for node ",n.name)
    return sequences

def run_ancestral(T, aln, reference_sequence=None, is_vcf=False, full_sequences=False, fill_overhangs=False,
                  infer_ambiguous=False, marginal=False, alphabet='nuc'):
    """
    ancestral nucleotide reconstruction using TreeTime
    """

    tt = ancestral_sequence_inference(tree=T, aln=aln, ref=reference_sequence if is_vcf else None, marginal=marginal,
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
    if reference_sequence:
        root_seq = reference_sequence
    else:
        root_seq = tt.sequence(T.root, as_string=True)

    mask = create_mask(is_vcf, tt, reference_sequence, aln)
    mutations = collect_mutations(tt, mask, character_map, reference_sequence, infer_ambiguous)
    sequences = {}
    if full_sequences:
        sequences = collect_sequences(tt, mask, reference_sequence, infer_ambiguous)

    # Combine the mutations & sequences into a single dict which downstream code
    # expects
    nodes = defaultdict(dict)
    for n in tt.tree.find_clades():
        name = n.name
        if name in mutations:
            nodes[name]['muts'] = mutations[name]
        if name in sequences:
            nodes[name]['sequence'] = sequences[name]

    return {'tt': tt,
            'root_seq': root_seq,
            'mutations': {"nodes": nodes, "mask": mask}}


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("ancestral", help=__doc__)

    input_group = parser.add_argument_group(
        "inputs",
        "Tree and sequences to use for ancestral reconstruction"
    )
    input_group.add_argument('--tree', '-t', required=True, help="prebuilt Newick")
    input_group.add_argument('--alignment', '-a', help="alignment in FASTA or VCF format")
    input_group_ref = input_group.add_mutually_exclusive_group()
    input_group_ref.add_argument('--vcf-reference', type=str, metavar='FASTA',
                                 help='[VCF alignment only] file of the sequence the VCF was mapped to.'
                                      ' Differences between this sequence and the inferred root will be reported as mutations on the root branch.')
    input_group_ref.add_argument('--root-sequence', type=str,metavar='FASTA/GenBank',
                                 help='[FASTA alignment only] file of the sequence that is used as root for mutation calling.'
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
    amino_acid_options_group.add_argument('--genes', nargs='+', action='extend', help="genes to translate (list or file containing list)")
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

    general_group = parser.add_argument_group(
        "general",
    )
    general_group.add_argument('--validation-mode', type=ValidationMode, choices=[mode for mode in ValidationMode], default=ValidationMode.ERROR,
                               help=validation_mode_help_message)

    return parser

def validate_arguments(args, is_vcf):
    """
    Check that provided arguments are compatible.
    Where possible we use argparse built-ins, but they don't cover everything we want to check.
    This checking shouldn't be used by downstream code to assume arguments exist, however by checking for
    invalid combinations up-front we can exit quickly.
    """
    aa_arguments = (args.annotation, args.genes, args.translations)
    if any(aa_arguments) and not all(aa_arguments):
        raise AugurError("For amino acid sequence reconstruction, you must provide an annotation file, a list of genes, and a template path to amino acid sequences.")

    if args.output_sequences and args.output_vcf:
        raise AugurError("Both sequence (fasta) and VCF output have been requested, but these are incompatible.")

    if is_vcf and args.output_sequences:
        raise AugurError("Sequence (fasta) output has been requested but the input alignment is VCF.")

    if not is_vcf and args.output_vcf:
        raise AugurError("VCF output has been requested but the input alignment is not VCF.")


def run(args):
    # check alignment type, set flags, read in if VCF
    is_vcf = is_filename_vcf(args.alignment)
    ref = None
    validate_arguments(args, is_vcf)

    try:
        T = read_tree(args.tree)
    except FileNotFoundError:
        raise AugurError(f"The provided tree file {args.tree!r} doesn't exist")
    except InvalidTreeError as error:
        raise AugurError(error)
    # Note that a number of other errors may be thrown by `read_tree` such as Bio.Phylo.NewickIO.NewickError

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
            raise AugurError("a reference Fasta is required with VCF-format alignments")
        compress_seq = read_vcf(args.alignment, args.vcf_reference)
        aln = compress_seq['sequences']
        ref = compress_seq['reference']
        vcf_metadata = compress_seq['metadata']
    else:
        aln = args.alignment
        ref = None
        if args.root_sequence:
            for fmt in ['fasta', 'genbank']:
                try:
                    ref = str(SeqIO.read(args.root_sequence, fmt).seq).upper()
                    break
                except:
                    pass
            if ref is None:
                raise AugurError(f"could not read root sequence from {args.root_sequence}")

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
    nuc_result = run_ancestral(T, aln, reference_sequence=ref if ref else None, is_vcf=is_vcf, fill_overhangs=not args.keep_overhangs,
                               full_sequences=full_sequences, marginal=args.inference, infer_ambiguous=infer_ambiguous, alphabet='nuc')
    anc_seqs = nuc_result['mutations']
    anc_seqs['reference'] = {'nuc': nuc_result['root_seq']}

    if anc_seqs.get("mask") is not None:
        anc_seqs["mask"] = "".join(['1' if x else '0' for x in anc_seqs["mask"]])

    anc_seqs['annotations'] = {'nuc': {'start': 1, 'end': len(anc_seqs['reference']['nuc']),
                                       'strand': '+', 'type': 'source'}}

    if not is_vcf and args.genes:
        genes = parse_genes_argument(args.genes)

        from .utils import load_features
        ## load features; only requested features if genes given
        features = load_features(args.annotation, genes)
        # Ensure the already-created nuc annotation coordinates match those parsed from the reference file
        if (features['nuc'].location.start+1 != anc_seqs['annotations']['nuc']['start'] or
            features['nuc'].location.end != anc_seqs['annotations']['nuc']['end']):
            raise AugurError(f"The 'nuc' annotation coordinates parsed from {args.annotation!r} ({features['nuc'].location.start+1}..{features['nuc'].location.end})"
                f" don't match the provided sequence data coordinates ({anc_seqs['annotations']['nuc']['start']}..{anc_seqs['annotations']['nuc']['end']}).")
        
        print("Read in {} features from reference sequence file".format(len(features)))
        for gene in genes:
            print(f"Processing gene: {gene}")
            fname = args.translations.replace("%GENE", gene)
            feat = features[gene]
            reference_sequence = str(feat.extract(Seq(ref)).translate()) if ref else None

            aa_result = run_ancestral(T, fname, reference_sequence=reference_sequence, is_vcf=is_vcf, fill_overhangs=not args.keep_overhangs,
                                        marginal=args.inference, infer_ambiguous=infer_ambiguous, alphabet='aa')
            if aa_result['tt'].data.full_length*3 != len(feat):
                raise AugurError(f"length of translated alignment for {gene} does not match length of reference feature."
                       " Please make sure that the annotation matches the translations.")

            for key, node in anc_seqs['nodes'].items():
                if 'aa_muts' not in node: node['aa_muts'] = {}
                node['aa_muts'][gene] = aa_result['mutations']['nodes'][key]['muts']

                # Add amino acid sequences to the root node of the tree.
                if key == T.root.name:
                    if "aa_sequences" not in node:
                        node["aa_sequences"] = {}

                    node["aa_sequences"][gene] = aa_result['tt'].sequence(T.root, as_string=True, reconstructed=True)

            anc_seqs['reference'][gene] = aa_result['root_seq']
            anc_seqs['annotations'].update(genome_features_to_auspice_annotation({gene: feat}, args.annotation))

            # Save ancestral amino acid sequences to FASTA.
            if args.output_translations:
                with open_file(args.output_translations.replace("%GENE", gene), "w") as oh:
                    for node in aa_result["tt"].tree.find_clades():
                        oh.write(f">{node.name}\n{aa_result['tt'].sequence(node, as_string=True, reconstructed=True)}\n")

    out_name = get_json_name(args, '.'.join(args.alignment.split('.')[:-1]) + '_mutations.json')
    # use NodeDataObject to perform validation on the file before it's written
    NodeDataObject(anc_seqs, out_name, args.validation_mode)

    write_json(anc_seqs, out_name)
    print("ancestral mutations written to", out_name, file=sys.stdout)

    if args.output_sequences:
        assert not is_vcf
        records = [
            SeqRecord(Seq(node_data["sequence"]), id=node_name, description="")
            for node_name, node_data in anc_seqs["nodes"].items()
        ]
        SeqIO.write(records, args.output_sequences, "fasta")
        print("ancestral sequences FASTA written to", args.output_sequences, file=sys.stdout)

    # output VCF including new ancestral seqs
    if args.output_vcf:
        assert is_vcf
        tree_dict = nuc_result['tt'].get_tree_dict(keep_var_ambigs=not infer_ambiguous)
        tree_dict['metadata'] = vcf_metadata
        write_vcf(tree_dict, args.output_vcf, anc_seqs['mask'])
        print("Mutations, including for ancestral nodes, exported as VCF to", args.output_vcf, file=sys.stdout)

    return 0
