"""
Infer ancestral sequences based on a tree.
"""

import os, shutil, time, json, sys
from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .utils import read_tree, InvalidTreeError, write_json
from treetime.vcf_utils import read_vcf, write_vcf
from collections import defaultdict

def ancestral_sequence_inference(tree=None, aln=None, ref=None, infer_gtr=True,
                                 marginal=False, fill_overhangs=True):
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

    Returns
    -------
    TreeAnc
        treetime.TreeAnc instance
    """

    from treetime import TreeAnc
    tt = TreeAnc(tree=tree, aln=aln, ref=ref, gtr='JC69',
                 fill_overhangs=fill_overhangs, verbose=1)

    # convert marginal (from args.inference) from 'joint' or 'marginal' to True or False
    bool_marginal = (marginal == "marginal")

    # only infer ancestral sequences, leave branch length untouched
    tt.infer_ancestral_sequences(infer_gtr=infer_gtr, marginal=bool_marginal)

    print("\nInferred ancestral sequence states using TreeTime:"
          "\n\tSagulenko et al. TreeTime: Maximum-likelihood phylodynamic analysis"
          "\n\tVirus Evolution, vol 4, https://academic.oup.com/ve/article/4/1/vex042/4794731\n")

    return tt

def collect_sequences_and_mutations(T, is_vcf=False):
    """iterates of the tree and produces dictionaries with
    mutations and sequences for each node.

    Parameters
    ----------
    T : Bio.Phylo.Tree
        Phylogenetic tree decorated with sequences and mutations as output by treetime.
    is_vcf : bool, optional
        specifies whether input alignment was vcf type (implying long genomes)

    Returns
    -------
    dict
        dictionary of mutations and sequences
    """
    data = defaultdict(dict)
    inc = 1 # convert python numbering to start-at-1
    for n in T.find_clades():
        if hasattr(n, "mutations"):
            mutations_attr = n.__getattribute__("mutations")
            data[n.name]['muts'] = [str(a)+str(int(pos)+inc)+str(d)
                                    for a,pos,d in mutations_attr]
    if not is_vcf:
        for n in T.find_clades():
            if hasattr(n, "sequence"):
                sequence_attr = n.__getattribute__("sequence")
                data[n.name]['sequence'] = ''.join(sequence_attr)
            else:
                data[n.name]['sequence'] = ''

    return data


def register_arguments(parser):
    parser.add_argument('--tree', '-t', required=True, help="prebuilt Newick")
    parser.add_argument('--alignment', '-a', help="alignment in fasta or VCF format")
    parser.add_argument('--output', '-o', type=str, help='name of JSON file to save mutations and ancestral sequences to')
    parser.add_argument('--output-node-data', type=str, help='name of JSON file to save mutations and ancestral sequences to')
    parser.add_argument('--output-sequences', type=str, help='name of FASTA file to save ancestral sequences to (FASTA alignments only)')
    parser.add_argument('--inference', default='joint', choices=["joint", "marginal"],
                                    help="calculate joint or marginal maximum likelihood ancestral sequence states")
    parser.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the VCF was mapped to')
    parser.add_argument('--output-vcf', type=str, help='name of output VCF file which will include ancestral seqs')
    parser.add_argument('--keep-ambiguous', action="store_true", default=False,
                                help='do not infer nucleotides at ambiguous (N) sites on tip sequences (leave as N). Always true for VCF input.')
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

    # Only allow recovery of ambig sites for Fasta-input if TreeTime is version 0.5.6 or newer
    # Otherwise it returns nonsense.
    from distutils.version import StrictVersion
    import treetime
    if args.keep_ambiguous and not is_vcf and StrictVersion(treetime.version) < StrictVersion('0.5.6'):
        print("ERROR: Keeping ambiguous sites for Fasta-input requires TreeTime version 0.5.6 or newer."+
                "\nYour version is "+treetime.version+
                "\nUpdate TreeTime or run without the --keep-ambiguous flag.")
        return 1

    tt = ancestral_sequence_inference(tree=T, aln=aln, ref=ref, marginal=args.inference,
                                      fill_overhangs = not(args.keep_overhangs))

    if is_vcf or args.keep_ambiguous:
        # TreeTime overwrites ambig sites on tips during ancestral reconst.
        # Put these back in tip sequences now, to avoid misleading
        tt.recover_var_ambigs()

    anc_seqs['nodes'] = collect_sequences_and_mutations(T, is_vcf)

    if args.output:
        anc_seqs_fname = args.output
        print("WARNING: the --output flag will be deprecated in the next major augur release. Use --output-node-data instead.", file=sys.stderr)
    elif args.output_node_data:
        anc_seqs_fname = args.output_node_data
    else:
        anc_seqs_fname = '.'.join(args.alignment.split('.')[:-1]) + '.anc_seqs.json'

    write_json(anc_seqs, anc_seqs_fname)
    print("ancestral mutations and sequences JSON written to",anc_seqs_fname, file=sys.stdout)

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
