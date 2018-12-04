"""
Annotate a tree with actual and inferred titer measurements.
"""

import json, os, sys
import numpy as np
from collections import defaultdict
from Bio import Phylo
from .utils import read_metadata, read_node_data, write_json


def register_arguments(parser):
    subparsers = parser.add_subparsers()
    sub_model = subparsers.add_parser('sub', help='substitution model')
    tree_model = subparsers.add_parser('tree', help='tree model')


    tree_model.add_argument('--titers', required=True, type=str, help="file with titer measurements")
    tree_model.add_argument('--tree', '-t', type=str, required=True, help="tree to perform fit titer model to")
    tree_model.add_argument('--output', '-o', type=str, help='JSON file to save titer model')
    tree_model.set_defaults(model = 'tree')

    sub_model.add_argument('--titers', required=True, type=str, help="file with titer measurements")
    sub_model.add_argument('--alignment', nargs='+', type=str, help="sequence to be used in the substitution model, supplied as fasta files")
    sub_model.add_argument('--gene-names', nargs='+', type=str, help="names of the sequences in the alignment, same order assumed")
    sub_model.add_argument('--output', '-o', type=str, help='JSON file to save titer model')
    sub_model.set_defaults(model = 'sub')


def load_alignments(sequence_files, gene_names):
    from Bio import AlignIO
    alignments = {}
    for fname, gene in zip(sequence_files, gene_names):
        alignments[gene] = AlignIO.read(fname, 'fasta')
    return alignments


def run(args):
    if args.model=='sub':
        infer_substitution_model(args)
    elif args.model=='tree':
        infer_tree_model(args)


def infer_substitution_model(args):
    from .titer_model import SubstitutionModel
    if not args.alignment:
        print('ERROR: substitution model requires an alignment. Please specify via --alignment')
        sys.exit(1)

    alignments = load_alignments(args.alignment, args.gene_names)

    TM_subs = SubstitutionModel(alignments, args.titers)
    TM_subs.prepare()
    TM_subs.train()

    # export the substitution model
    subs_model = {'titers':TM_subs.compile_titers(),
                  'potency':TM_subs.compile_potencies(),
                  'avidity':TM_subs.compile_virus_effects(),
                  'substitution':TM_subs.compile_substitution_effects()}
    write_json(subs_model, args.output)

    print("\nInferred titer model of type 'SubstitutionModel' using augur:"
          "\n\tNeher et al. Prediction, dynamics, and visualization of antigenic phenotypes of seasonal influenza viruses."
          "\n\tPNAS, vol 113, 10.1073/pnas.1525578113\n")
    print("results written to", args.output)


def infer_tree_model(args):
    T = Phylo.read(args.tree, 'newick')
    from .titer_model import TreeModel
    TM_tree = TreeModel(T, args.titers)
    TM_tree.prepare()
    TM_tree.train()

    # export the tree model
    tree_model = {'titers':TM_tree.compile_titers(),
                  'potency':TM_tree.compile_potencies(),
                  'avidity':TM_tree.compile_virus_effects(),
                  'nodes':{n.name:{"dTiter": n.dTiter, "cTiter":n.cTiter}
                              for n in T.find_clades()}}
    write_json(tree_model, args.output)
    print("\nInferred titer model of type 'TreeModel' using augur:"
          "\n\tNeher et al. Prediction, dynamics, and visualization of antigenic phenotypes of seasonal influenza viruses."
          "\n\tPNAS, vol 113, 10.1073/pnas.1525578113\n")
    print("results written to", args.output)
