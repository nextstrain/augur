import numpy as np
from collections import defaultdict
import json, os
from Bio import Phylo
from .utils import read_metadata, read_node_data, write_json

def run(args):
    T = Phylo.read(args.tree, 'newick')

    TM_tree, TM_subs = None, None
    if args.titer_model == "tree":
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

    if args.titer_model == "substitution":
        from .titer_model import SubstitutionModel
        TM_subs = SubstitutionModel(T, args.titers)
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
