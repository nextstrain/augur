"""
Infer ancestral sequences based on a tree.
"""

import os, shutil, time, json, sys
from Bio import Phylo
from .utils import write_json
from treetime.vcf_utils import read_vcf, write_vcf
from collections import defaultdict

def ancestral_sequence_inference(tree=None, aln=None, ref=None, infer_gtr=True,
                                 marginal=False):
    from treetime import TreeAnc
    tt = TreeAnc(tree=tree, aln=aln, ref=ref, gtr='JC69', verbose=1)

    # convert marginal (from args.inference) from 'joint' or 'marginal' to True or False
    bool_marginal = (marginal == "marginal")

    # only infer ancestral sequences, leave branch length untouched
    tt.infer_ancestral_sequences(infer_gtr=infer_gtr, marginal=bool_marginal)

    print("\nInferred ancestral sequence states using TreeTime:"
          "\n\tSagulenko et al. TreeTime: Maximum-likelihood phylodynamic analysis"
          "\n\tVirus Evolution, vol 4, https://academic.oup.com/ve/article/4/1/vex042/4794731\n")

    return tt

def collect_sequences_and_mutations(T, is_vcf=False):
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
    parser.add_argument('--output', '-o', type=str, help='file name to save mutations and ancestral sequences to')
    parser.add_argument('--inference', default='joint', choices=["joint", "marginal"],
                                    help="calculate joint or marginal maximum likelihood ancestral sequence states")
    parser.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the VCF was mapped to')
    parser.add_argument('--output-vcf', type=str, help='name of output VCF file which will include ancestral seqs')


def run(args):
    # check alignment type, set flags, read in if VCF
    is_vcf = False
    ref = None
    anc_seqs = {}
    # check if tree is provided and can be read
    for fmt in ["newick", "nexus"]:
        try:
            T = Phylo.read(args.tree, fmt)
            break
        except:
            pass
    if T is None:
        print("ERROR: reading tree from %s failed."%args.tree)
        return 1

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

    tt = ancestral_sequence_inference(tree=T, aln=aln, ref=ref, marginal=args.inference)

    if is_vcf:
        # TreeTime overwrites ambig sites on tips during ancestral reconst.
        # Put these back in tip sequences now, to avoid misleading
        tt.recover_var_ambigs()

    anc_seqs['nodes'] = collect_sequences_and_mutations(T, is_vcf)

    if args.output:
        anc_seqs_fname = args.output
    else:
        anc_seqs_fname = '.'.join(args.alignment.split('.')[:-1]) + '.anc_seqs.json'

    anc_seqs_success = write_json(anc_seqs, anc_seqs_fname)
    print("ancestral sequences written to",anc_seqs_fname, file=sys.stdout)

    # If VCF, output VCF including new ancestral seqs
    if is_vcf:
        if args.output_vcf:
            vcf_fname = args.output_vcf
        else:
            vcf_fname = '.'.join(args.alignment.split('.')[:-1]) + '.vcf'
        write_vcf(tt.get_tree_dict(keep_var_ambigs=True), vcf_fname)
        print("ancestral sequences as vcf-file written to",vcf_fname, file=sys.stdout)

    if anc_seqs_success:
        return 0
    else:
        return 1
