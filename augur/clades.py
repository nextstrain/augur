import os, sys
import numpy as np
from Bio import SeqIO, SeqFeature, Seq, SeqRecord, Phylo
from .utils import read_node_data, write_json
from collections import defaultdict

def read_in_clade_definitions(clade_file):
    '''
    Reads in tab-seperated file that defines clades by amino-acid

    Format:
    clade	gene	site	aa
    Clade_2	embC	940	S
    Clade_2	Rv3463	192	K
    Clade_3	Rv2209	432	I
    '''
    import pandas as pd

    clades = defaultdict(lambda:defaultdict(list))

    df = pd.read_csv(clade_file, sep='\t' if clade_file.endswith('.tsv') else ',')
    for mi, m in df.iterrows():
        clades[m.clade][m.gene].append((m.site,m.aa))

    return clades

def assign_clades(clade_designations, aa_muts, tree):
    clades = {}
    isClade = False
    for n in tree.get_nonterminals():
        n_muts = aa_muts[n.name]['aa_muts']
        if n.name not in clades:
            clades[n.name]={"clade_membership": "Unassigned"}
        for clade, definition in clade_designations.items():
            isClade = False
            for gene, muts in definition.items():
                if (gene in n_muts and n_muts[gene] != []):
                    prs_muts = [(int(tu[1:-1]), tu[-1]) for tu in n_muts[gene]] #get mutations in right format
                    if all([mut in prs_muts for mut in muts]):
                        isClade = True
                    else:
                        isClade = False
                        break
                else:
                    isClade = False
                    break
            if isClade:
                clades[n.name] = {"clade_annotation":clade, "clade_membership": clade}

        for c in n:
            clades[c.name]={"clade_membership": clades[n.name]["clade_membership"] }

    return clades

def run(args):
    ## read tree and data, if reading data fails, return with error code
    tree = Phylo.read(args.tree, 'newick')
    node_data = read_node_data(args.amino_acids, args.tree)
    if node_data is None:
        print("ERROR: could not read node data (incl sequences)")
        return -1

    clade_designations = read_in_clade_definitions(args.clades)

    aa_muts = node_data['nodes']

    clades = assign_clades(clade_designations, aa_muts, tree)

    write_json({'nodes':clades}, args.output)
    print("clades written to", args.output, file=sys.stdout)
