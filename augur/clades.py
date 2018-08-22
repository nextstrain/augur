import os, sys
import numpy as np
from Bio import SeqIO, SeqFeature, Seq, SeqRecord, Phylo
from .utils import read_node_data, write_json
from collections import defaultdict

def read_in_clade_definitions(clade_file):
    '''
    Reads in tab-seperated file that defines clades by amino-acid.

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


def is_node_in_clade(clade_mutations, node_muts):
    '''
    Determines whether a node contains all mutations that define a clade
    '''
    isClade = False
    for gene, muts in clade_mutations.items():
        if (gene in node_muts and node_muts[gene] != []):
            prs_muts = [(int(tu[1:-1]), tu[-1]) for tu in node_muts[gene]] #get mutations in right format
            if all([mut in prs_muts for mut in muts]):
                isClade = True
            else:
                return False
        else:
            return False

    return isClade


def assign_clades(clade_designations, muts, tree):
    '''
    Ensures all nodes have an entry (or auspice doesn't display nicely), tests each node
    to see if it's the first member of a clade (assigns 'clade_annotation'), and sets
    all nodes's clade_membership to the value of their parent. This will change if later found to be
    the first member of a clade.
    '''
    clades = {}
    for n in tree.get_nonterminals(order = 'preorder'):
        n_muts = {}
        if 'aa_muts' in muts[n.name]:
            n_muts = muts[n.name]['aa_muts']
        if 'muts' in muts[n.name]:
            n_muts['nuc'] = muts[n.name]['muts'] # Put nuc mutations in with 'nuc' as the 'gene' so all can be searched together

        if n.name not in clades: # This ensures every node gets an entry - otherwise auspice doesn't display nicely
            clades[n.name]={"clade_membership": "Unassigned"}
        for clade, definition in clade_designations.items():
            if is_node_in_clade(definition, n_muts):
                clades[n.name] = {"clade_annotation":clade, "clade_membership": clade}

        # Ensures each node is set to membership of their parent initially (unless changed later in tree traversal)
        for c in n:
            clades[c.name]={"clade_membership": clades[n.name]["clade_membership"] }

    return clades


def run(args):
    ## read tree and data, if reading data fails, return with error code
    tree = Phylo.read(args.tree, 'newick')
    node_data = read_node_data(args.mutations, args.tree)
    if node_data is None:
        print("ERROR: could not read node data (incl sequences)")
        return -1

    clade_designations = read_in_clade_definitions(args.clades)

    muts = node_data['nodes']

    clades = assign_clades(clade_designations, muts, tree)

    write_json({'nodes':clades}, args.output)
    print("clades written to", args.output, file=sys.stdout)
