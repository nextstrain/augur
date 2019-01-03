"""
Assign clades to nodes in a tree based on amino-acid or nucleotide signatures.
"""

import sys
from Bio import Phylo
import pandas as pd
from .utils import get_parent_name_by_child_name_for_tree, read_node_data, write_json

def read_in_clade_definitions(clade_file):
    '''
    Reads in tab-seperated file that defines clades by amino acid or nucleotide mutations

    Format:
    clade	gene	site alt
    Clade_1	ctpE    81  D
    Clade_2	nuc 30642   T
    Clade_3	nuc 444296  A
    Clade_4	pks8    634 T
    '''

    clades = {}

    df = pd.read_csv(clade_file, sep='\t' if clade_file.endswith('.tsv') else ',')
    for index, row in df.iterrows():
        allele = (row.gene, row.site, row.alt)
        if row.clade in clades:
            clades[row.clade].append(allele)
        else:
            clades[row.clade] = [allele]

    return clades

def get_node_alleles(node_name, muts, parents):
    '''
    Retrieve the set of mutations leading to a node
    This includes node-specific mutations as well as mutations
    along its full line of descent
    '''
    node_alleles = []
    sites_encountered = set()
    focal_node = node_name
    while focal_node is not None:
        if 'aa_muts' in muts[focal_node]:
            for gene in muts[focal_node]['aa_muts']:
                for mut in muts[focal_node]['aa_muts'][gene]:
                    site = (gene, int(mut[1:-1]))
                    if site not in sites_encountered:
                        sites_encountered.add(site)
                        allele = (gene, int(mut[1:-1]), mut[-1])
                        node_alleles.append(allele)
        if 'muts' in muts[focal_node]:
            for mut in muts[focal_node]['muts']:
                site = (gene, int(mut[1:-1]))
                if site not in sites_encountered:
                    sites_encountered.add(site)
                    allele = ('nuc', int(mut[1:-1]), mut[-1])
                    node_alleles.append(allele)

        focal_node = parents.get(focal_node)

    return node_alleles

def is_node_in_clade(clade_alleles, node_alleles):
    '''
    Determines whether a node contains all mutations that define a clade
    '''
    is_clade = False
    # mutations are stored on nodes in format 'R927H' matching the clade definitions
    # if all of the clade-defining mutations are in the node mutations, it's part of this clade
    if all([ clade_allele in node_alleles for clade_allele in clade_alleles ]):
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
    parents = get_parent_name_by_child_name_for_tree(tree)

    # first pass to set all nodes to unassigned as precaution to ensure attribute is set
    for node in tree.find_clades(order = 'preorder'):
        clade_membership[node.name] = {'clade_membership': 'unassigned'}

    # second pass to assign 'clade_annotation' to basal nodes within each clade
    # if multiple nodes match, assign annotation to largest
    # otherwise occasional unwanted cousin nodes get assigned the annotation
    for clade_name, clade_alleles in clade_designations.items():
        node_counts = {}
        for node in tree.find_clades(order = 'preorder'):
            node_alleles = get_node_alleles(node.name, all_muts, parents)
            if is_node_in_clade(clade_alleles, node_alleles):
                node_counts[node.name] = node.count_terminals()
        sorted_nodes = list(sorted(node_counts.items(), key=lambda x: x[1], reverse=True))
        if len(sorted_nodes) > 0:
            target_node = sorted_nodes[0][0]
            clade_membership[target_node] = {'clade_annotation': clade_name, 'clade_membership': clade_name}

    # third pass to propagate 'clade_membership'
    # don't propagate if encountering 'clade_annotation'
    for node in tree.find_clades(order = 'preorder'):
        for child in node:
            if 'clade_annotation' not in clade_membership[child.name]:
                clade_membership[child.name]['clade_membership'] = clade_membership[node.name]['clade_membership']

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
