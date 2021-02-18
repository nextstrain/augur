"""
Infer ancestral sequences based on a tree.
"""

import os, shutil, time, json, sys
import numpy as np
from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .utils import read_tree, InvalidTreeError, write_json, get_json_name
from treetime.vcf_utils import read_vcf, write_vcf
from collections import defaultdict

def ancestral_sequence_inference(tree=None, aln=None, ref=None, infer_gtr=True,
                                 marginal=False, fill_overhangs=True, infer_tips=False):
    """infer ancestral sequences using TreeTime

    Parameters
    ----------
    tree : Bio.Phylo tree or str
        tree or filename of tree
    aln : Bio.Align.MultipleSeqAlignment or str
        alignment or filename of alignment
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

    Returns
    -------
    TreeAnc
        treetime.TreeAnc instance
    """

    from treetime import TreeAnc, version as treetime_version
    print(f"augur ancestral is using TreeTime version {treetime_version}")

    tt = TreeAnc(tree=tree, aln=aln, ref=ref, gtr='JC69',
                 fill_overhangs=fill_overhangs, verbose=1)

    # convert marginal (from args.inference) from 'joint' or 'marginal' to True or False
    bool_marginal = (marginal == "marginal")

    # only infer ancestral sequences, leave branch length untouched
    tt.infer_ancestral_sequences(infer_gtr=infer_gtr, marginal=bool_marginal,
                                 reconstruct_tip_states=infer_tips)

    print("\nInferred ancestral sequence states using TreeTime:"
          "\n\tSagulenko et al. TreeTime: Maximum-likelihood phylodynamic analysis"
          "\n\tVirus Evolution, vol 4, https://academic.oup.com/ve/article/4/1/vex042/4794731\n")

    return tt

def collect_mutations_and_sequences(tt, infer_tips=False, full_sequences=False, character_map=None, mask_ambiguous=True):
    """iterates of the tree and produces dictionaries with
    mutations and sequences for each node.

    Parameters
    ----------
    tt : treetime
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

    mask=None
    if full_sequences:
        if mask_ambiguous:
            # Identify sites for which every terminal sequence is ambiguous.
            # These sites will be masked to prevent rounding errors in the
            # maximum likelihood inference from assigning an arbitrary
            # nucleotide to sites at internal nodes.
            ambiguous_count = np.zeros(tt.sequence_length, dtype=int)
            for n in tt.tree.get_terminals():
                ambiguous_count += np.array(tt.sequence(n,reconstructed=False, as_string=False)==tt.gtr.ambiguous, dtype=int)
            mask = ambiguous_count==len(tt.tree.get_terminals())
        else:
            mask = np.zeros(tt.sequence_length, dtype=bool)

        for n in tt.tree.find_clades():
            try:
                tmp = tt.sequence(n,reconstructed=infer_tips, as_string=False)
                tmp[mask] = tt.gtr.ambiguous
                data[n.name]['sequence'] = "".join(tmp)
            except:
                print("No sequence available for node ",n.name)

    return {"nodes": data, "mask": mask}


def register_arguments(parser):
    parser.add_argument('--tree', '-t', required=True, help="prebuilt Newick")
    parser.add_argument('--alignment', '-a', help="alignment in fasta or VCF format")
    parser.add_argument('--output-node-data', type=str, help='name of JSON file to save mutations and ancestral sequences to')
    parser.add_argument('--output-sequences', type=str, help='name of FASTA file to save ancestral sequences to (FASTA alignments only)')
    parser.add_argument('--inference', default='joint', choices=["joint", "marginal"],
                                    help="calculate joint or marginal maximum likelihood ancestral sequence states")
    parser.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the VCF was mapped to')
    parser.add_argument('--output-vcf', type=str, help='name of output VCF file which will include ancestral seqs')
    ambiguous = parser.add_mutually_exclusive_group()
    ambiguous.add_argument('--keep-ambiguous', action="store_true",
                                help='do not infer nucleotides at ambiguous (N) sites on tip sequences (leave as N).')
    ambiguous.add_argument('--infer-ambiguous', action="store_true", default=True,
                                help='infer nucleotides at ambiguous (N,W,R,..) sites on tip sequences and replace with most likely state.')
    parser.add_argument('--keep-overhangs', action="store_true", default=False,
                                help='do not infer nucleotides for gaps (-) on either side of the alignment')

def run(args):
    # check alignment type, set flags, read in if VCF
    is_vcf = False
    ref = None
    anc_seqs = {}

    try:
        T = read_tree(args.tree)
    except (FileNotFoundError, InvalidTreeError) as error:
        print("ERROR: %s" % error, file=sys.stderr)
        return 1

    import numpy as np
    missing_internal_node_names = [n.name is None for n in T.get_nonterminals()]
    if np.all(missing_internal_node_names):
        print("\n*** WARNING: Tree has no internal node names!")
        print("*** Without internal node names, ancestral sequences can't be linked up to the correct node later.")
        print("*** If you want to use 'augur export' or `augur translate` later, re-run this command with the output of 'augur refine'.")
        print("*** If you haven't run 'augur refine', you can add node names to your tree by running:")
        print("*** augur refine --tree %s --output-tree <filename>.nwk"%(args.tree) )
        print("*** And use <filename>.nwk as the tree when running 'ancestral', 'translate', and 'traits'")

    if any([args.alignment.lower().endswith(x) for x in ['.vcf', '.vcf.gz']]):
        if not args.vcf_reference:
            print("ERROR: a reference Fasta is required with VCF-format alignments")
            return 1

        compress_seq = read_vcf(args.alignment, args.vcf_reference)
        aln = compress_seq['sequences']
        ref = compress_seq['reference']
        is_vcf = True
    else:
        aln = args.alignment

    # Enfore treetime 0.7 or later
    from distutils.version import StrictVersion
    import treetime
    if StrictVersion(treetime.version) < StrictVersion('0.7.0'):
        print("ERROR: this version of augur requires TreeTime 0.7 or later.")
        return 1

    # Infer ambiguous bases if the user has requested that we infer them (either
    # explicitly or by default) and the user has not explicitly requested that
    # we keep them.
    infer_ambiguous = args.infer_ambiguous and not args.keep_ambiguous

    tt = ancestral_sequence_inference(tree=T, aln=aln, ref=ref, marginal=args.inference,
                                      fill_overhangs = not(args.keep_overhangs),
                                      infer_tips = infer_ambiguous)

    character_map = {}
    for x in tt.gtr.profile_map:
        if tt.gtr.profile_map[x].sum()==tt.gtr.n_states:
            # TreeTime treats all characters that are not valid IUPAC nucleotide chars as fully ambiguous
            # To clean up auspice output, we map all those to 'N'
            character_map[x] = 'N'
        else:
            character_map[x] = x

    anc_seqs.update(collect_mutations_and_sequences(tt, full_sequences=not is_vcf,
                          infer_tips=infer_ambiguous, character_map=character_map))
    # add reference sequence to json structure. This is the sequence with
    # respect to which mutations on the tree are defined.
    if is_vcf:
        anc_seqs['reference'] = {"nuc":compress_seq['reference']}
    else:
        root_seq = tt.sequence(T.root, as_string=False)
        if anc_seqs.get("mask") is not None:
            root_seq[anc_seqs['mask']] = tt.gtr.ambiguous
        anc_seqs['reference'] = {"nuc": ''.join(root_seq)}

    if anc_seqs.get("mask") is not None:
        anc_seqs["mask"] = "".join(['1' if x else '0' for x in anc_seqs["mask"]])

    out_name = get_json_name(args, '.'.join(args.alignment.split('.')[:-1]) + '_mutations.json')
    write_json(anc_seqs, out_name)
    print("ancestral mutations written to", out_name, file=sys.stdout)

    if args.output_sequences:
        if args.output_vcf:
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
        write_vcf(tt.get_tree_dict(keep_var_ambigs=True), vcf_fname)
        print("ancestral sequences as vcf-file written to",vcf_fname, file=sys.stdout)

    return 0
