"""
Assign clades to nodes in a tree based on amino-acid or nucleotide signatures.
"""

import os, sys
import numpy as np
from Bio import SeqIO, SeqFeature, Seq, SeqRecord, Phylo
from .utils import read_node_data, write_json
from collections import defaultdict

def read_in_clade_definitions(clade_file):
    '''
    Reads in tab-seperated file that defines clades by amino-acid or nucleotide

    Format:
    clade	gene	site	alt
    Clade_2	embC	940	S
    Clade_2	Rv3463	192	K
    Clade_3	Rv2209	432	I
    Clade_4	nuc	27979	G
    '''
    import pandas as pd

    clades = defaultdict(lambda:defaultdict(list))

    df = pd.read_csv(clade_file, sep='\t' if clade_file.endswith('.tsv') else ',')
    for mi, m in df.iterrows():
        clades[m.clade][m.gene].append((m.site,m.alt))

    return clades


def is_node_in_clade(clade_mutations, node_muts):
    '''
    Determines whether a node contains all mutations that define a clade
    '''
    is_clade = False
    # cycle through each gene in the clade definition and see if this node has the specified mutations
    for gene, defining_mutations in clade_mutations.items():
        #check the node has the gene in question and that it contains any mutations
        if (gene in node_muts and node_muts[gene]):
            # mutations are stored on nodes in format 'R927H', but in clade defining_mutations as (927, 'H')
            # So convert node mutations ('mutation') to format (927, 'H') here, for easy comparison
            formatted_mutations = [(int(mutation[1:-1]), mutation[-1]) for mutation in node_muts[gene]]
            # If all of the clade-defining mutations are in the node mutations, it's part of this clade.
            if all([mut in formatted_mutations for mut in defining_mutations]):
                is_clade = True
            else:
                return False
        else:
            return False

    return is_clade


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
            clades[n.name]={"clade_membership": "unassigned"}
        for clade, definition in clade_designations.items():
            if is_node_in_clade(definition, n_muts):
                clades[n.name] = {"clade_annotation":clade, "clade_membership": clade}

        # Ensures each node is set to membership of their parent initially (unless changed later in tree traversal)
        for c in n:
            clades[c.name]={"clade_membership": clades[n.name]["clade_membership"] }

    return clades


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

    clade_designations = read_in_clade_definitions(args.clades)

    muts = node_data['nodes']

    clades = assign_clades(clade_designations, muts, tree)

    write_json({'nodes':clades}, args.output)
    print("clades written to", args.output, file=sys.stdout)
