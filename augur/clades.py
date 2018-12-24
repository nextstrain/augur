"""
Assign clades to nodes in a tree based on amino-acid or nucleotide signatures.
"""

import sys
from Bio import Phylo
import pandas as pd
from .utils import read_node_data, write_json

def read_in_clade_definitions(clade_file):
    '''
    Reads in tab-seperated file that defines clades by amino acid or nucleotide mutations

    Format:
    clade	gene	mutation
    Clade_1	ctpE	G81D
    Clade_2	nuc	C30642T
    Clade_3	nuc	C444296A
    Clade_4	pks8	A634T
    '''

    clades = {}

    df = pd.read_csv(clade_file, sep='\t' if clade_file.endswith('.tsv') else ',')
    for index, row in df.iterrows():
        pair = (row.gene, row.mutation)
        if row.clade in clades:
            clades[row.clade].append(pair)
        else:
            clades[row.clade] = [pair]
    return clades

def all_parents(tree):
    '''
    Return dictionary mapping child node names to parent node names
    '''
    parents = {}
    for clade in tree.find_clades(order='level'):
        for child in clade:
            parents[child.name] = clade.name
    return parents

def get_node_mutations(node_name, muts, parents):
    '''
    Retrieve the set of mutations leading to a node
    This includes node-specific mutations as well as mutations
    along its full line of descent
    '''
    node_muts = []
    focal_node = node_name
    while focal_node is not None:
        if 'aa_muts' in muts[focal_node]:
            for gene in muts[focal_node]['aa_muts']:
                for mut in muts[focal_node]['aa_muts'][gene]:
                    pair = (gene, mut)
                    node_muts.append(pair)
        if 'muts' in muts[focal_node]:
            for mut in muts[focal_node]['muts']:
                pair = ('nuc', mut)
                node_muts.append(pair)
        focal_node = parents.get(focal_node, None)
    return node_muts

def is_node_in_clade(clade_muts, node_muts):
    '''
    Determines whether a node contains all mutations that define a clade
    '''
    is_clade = False
    # mutations are stored on nodes in format 'R927H' matching the clade definitions
    # if all of the clade-defining mutations are in the node mutations, it's part of this clade
    if all([ clade_mut in node_muts for clade_mut in clade_muts ]):
        is_clade = True
    return is_clade

def assign_clades(clade_designations, all_muts, tree):
    '''
    Ensures all nodes have an entry (or auspice doesn't display nicely), tests each node
    to see if it's the first member of a clade (assigns 'clade_annotation'), and sets
    all nodes's clade_membership to the value of their parent. This will change if later found to be
    the first member of a clade.
    '''
    clade_membership = {}
    # first pass to set all nodes to unassigned
    for node in tree.get_nonterminals(order = 'preorder'):
        clade_membership[node.name] = {"clade_membership": "unassigned"}

    parents = all_parents(tree)

    for clade_name, clade_muts in clade_designations.items():
        first_instance_of_clade = True
        for node in tree.get_nonterminals(order = 'preorder'):
            node_muts = get_node_mutations(node.name, all_muts, parents)
            if is_node_in_clade(clade_muts, node_muts):
                if first_instance_of_clade:
                    clade_membership[node.name] = {"clade_annotation": clade_name, "clade_membership": clade_name}
                    first_instance_of_clade = False
                else:
                    clade_membership[node.name] = {"clade_membership": clade_name}
                # Ensures each node is set to membership of their parent initially (unless changed later in tree traversal)
                for c in node:
                    clade_membership[c.name] = {"clade_membership": clade_name}

    return clade_membership


def register_arguments(parser):
    parser.add_argument('--tree', help="prebuilt Newick -- no tree will be built if provided")
    parser.add_argument('--mutations', nargs='+', help='JSON(s) containing ancestral and tip nucleotide and/or amino-acid mutations ')
    parser.add_argument('--clades', type=str, help='TSV file containing clade definitions by amino-acid')
    parser.add_argument('--output', type=str, help="name of JSON files for clades")


def run(args):
    ## read tree and data, if reading data fails, return with error code
    tree = Phylo.read(args.tree, 'newick')
    node_data = read_node_data(args.mutations, args.tree)
    if node_data is None:
        print("ERROR: could not read node data (incl sequences)")
        return 1
    all_muts = node_data['nodes']

    clade_designations = read_in_clade_definitions(args.clades)

    clade_membership = assign_clades(clade_designations, all_muts, tree)

    write_json({'nodes': clade_membership}, args.output)
    print("clades written to", args.output, file=sys.stdout)
