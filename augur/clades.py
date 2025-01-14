"""
Assign clades to nodes in a tree based on amino-acid or nucleotide signatures.

Nodes which are members of a clade are stored via
<OUTPUT_NODE_DATA> → nodes → <node_name> → clade_membership
and if this file is used in `augur export v2` these will automatically become a coloring.

The basal nodes of each clade are also given a branch label which is stored via
<OUTPUT_NODE_DATA> → branches → <node_name> → labels → clade.

The keys "clade_membership" and "clade" are customisable via command line arguments.
"""

import sys
from Bio import Phylo
import pandas as pd
import numpy as np
from collections import defaultdict
import networkx as nx
from itertools import islice
from .argparse_ import ExtendOverwriteDefault
from .errors import AugurError
from .io.file import PANDAS_READ_CSV_OPTIONS
from argparse import SUPPRESS
from .utils import get_parent_name_by_child_name_for_tree, read_node_data, write_json, get_json_name
from .argparse_ import add_validation_arguments

UNASSIGNED = 'unassigned'

def read_in_clade_definitions(clade_file):
    '''
    Reads in tab-seperated file that defines clades by amino acid or nucleotide mutations

    Inheritance is allowed, but needs to be acyclic. Alleles can be overwritten by inheriting clades.

    Sites are 1 indexed in the file, and are converted to 0 indexed in the output

    Empty lines are ignored, comments after # are ignored

    Format::

        clade      gene    site     alt
        Clade_1    ctpE    81       D
        Clade_2    nuc     30642    T
        Clade_3    nuc     444296   A
        Clade_3    S       1        P
        # Clade_4 inherits from Clade_3
        Clade_4    clade   Clade_3
        Clade_4    pks8    634      T
        # Inherited allele can be overwritten
        Clade_4    S       1        L

    Parameters
    ----------
    clade_file : str
        meta data file

    Returns
    -------
    dict
        clade definitions as :code:`{clade_name:[(gene, site, allele),...]}`
    '''

    clades = defaultdict(lambda: defaultdict(str))
    df = pd.read_csv(
        clade_file,
        sep='\t' if clade_file.endswith('.tsv') else ',',
        comment='#',
        na_filter=False,
        **PANDAS_READ_CSV_OPTIONS,
    )

    clade_inheritance_rows = df[df['gene'] == 'clade']

    # Identify clades that inherit more than once
    clades_with_multiple_inheritance = clade_inheritance_rows[clade_inheritance_rows.duplicated(subset=["clade"])]['clade'].tolist()
    if len(clades_with_multiple_inheritance) > 0:
        raise ValueError(f"Clades {clades_with_multiple_inheritance} have multiple inheritance, that's not allowed")

    # Identify clades that inherit from non-existent clades
    missing_parent_clades = set(clade_inheritance_rows['site']) - set(df["clade"])
    if len(missing_parent_clades) > 0:
        raise ValueError(f"Clades {missing_parent_clades} are inherited from but are not defined")


    G = nx.DiGraph()

    # Use integer 0 as root so as not to conflict with any string clade names
    # String '0' can still be used this way
    root = 0

    # Skip rows that are missing a clade name.
    defined_clades = (clade for clade in df.clade.unique() if clade != '')

    # For every clade, add edge from root as default
    # This way all clades can be reached by traversal
    for clade in defined_clades:
        G.add_edge(root, clade)

    # Build inheritance graph
    # For clades that inherit, disconnect from root
    # Add edge from parent
    for _, row in clade_inheritance_rows.iterrows():
        G.remove_edge(root, row.clade)
        G.add_edge(row.site, row.clade)

    if not nx.is_directed_acyclic_graph(G):
        raise ValueError(f"Clade definitions contain cycles {list(nx.simple_cycles(G))}")

    # Traverse graph top down, so that children can inherit from parents and grandparents
    # Topological sort ensures parents are visited before children
    # islice is used to skip the root node (which has no parent)
    for clade in islice(nx.topological_sort(G),1,None):
        # Get name of parent clade
        # G.predecessors(clade) returns iterator, thus next() necessary
        # despite the fact that there should only be one parent
        parent_clade = next(G.predecessors(clade))
        # Inheritance from parents happens here
        # Allele dict is initialized with alleles from parent
        clades[clade] = clades[parent_clade].copy()
        for _, row in df[(df.clade == clade) & (df.gene != 'clade')].iterrows():
            # Overwrite of parent alleles is possible and happens here
            clades[clade][(row.gene, int(row.site)-1)] = row.alt

    # Convert items from dict[str, dict[(str,int),str]] to dict[str, list[(str,int,str)]]
    clades = {
        clade: [
            gene_site + (alt,)
            for gene_site, alt in clade_definition.items()
        ]
        for clade, clade_definition in clades.items()
        # If clause avoids root (helper) from being emmitted
        if clade != root
    }

    if not len(clades.keys()):
        raise AugurError(f"No clades were defined in {clade_file}")

    return clades


def is_node_in_clade(clade_alleles, node, root_sequence):
    '''
    Determines whether a node matches the clade definition based on sequence
    For any condition, will first look in mutations stored in node.sequences,
    then check whether a reference sequence is available, and other reports 'non-match'

    Parameters
    ----------
    clade_alleles : list
        list of clade defining alleles (typically supplied from the input TSV)
    node : Bio.Phylo.BaseTree.Clade
        node to check, assuming sequences (as mutations) are attached to node
        node.sequences specifies nucleotides/codons which are newly observed on this node
        i.e. they are the result of a mutation observed on the branch leading to this node
    root_sequence : dict
        {geneName: observed root sequence (list)}

    Returns
    -------
    bool
        True if in clade

    '''
    conditions = []
    for gene, pos, clade_state in clade_alleles:
        if gene in node.sequences and pos in node.sequences[gene]:
            state = node.sequences[gene][pos]
        elif root_sequence and gene in root_sequence:
            try:
                state = root_sequence[gene][pos]
            except IndexError:
                raise AugurError(f"A clade definition specifies {{{gene},{pos+1},{clade_state}}} which \
is beyond the bounds of the supplied root sequence for {gene} (length {len(root_sequence[gene])})")
        else:
            state = ''

        conditions.append(state==clade_state)

    return all(conditions)

def ensure_no_multiple_mutations(all_muts):
    multiples = []

    for name,node in all_muts.items():
        nt_positions = [int(mut[1:-1])-1  for mut in node.get('muts', [])]
        if len(set(nt_positions))!=len(nt_positions):
            multiples.append(f"Node {name} (nuc)")
        for gene in node.get('aa_muts', {}):
            aa_positions = [int(mut[1:-1])-1 for mut in node['aa_muts'][gene]]
            if len(set(aa_positions))!=len(aa_positions):
                multiples.append(f"Node {name} ({gene})")

    if multiples:
        raise AugurError(f"Multiple mutations at the same position on a single branch were found: {', '.join(multiples)}")

def assign_clades(clade_designations, all_muts, tree, ref=None):
    '''
    Ensures all nodes have an entry (or auspice doesn't display nicely), tests each node
    to see if it's the first member of a clade (this is the label), and sets the membership of each
    node to the value of their parent. This will change if later found to be
    the first member of a clade.

    Parameters
    ----------
    clade_designations :     dict
        clade definitions as :code:`{clade_name:[(gene, site, allele),...]}`
    all_muts : dict
        mutations in each node
    tree : Bio.Phylo.BaseTree.Tree
        phylogenetic tree to process
    ref : str or list, optional
        reference sequence to look up state when not mutated

    Returns
    -------
    (dict, dict)
        [0]: mapping of node to clade membership (where applicable)
        [1]: mapping of node to clade label (where applicable)
    '''

    clade_membership = {}
    clade_labels = {}
    parents = get_parent_name_by_child_name_for_tree(tree)

    # first pass to set all nodes to unassigned as precaution to ensure attribute is set
    for node in tree.find_clades(order = 'preorder'):
        clade_membership[node.name] = UNASSIGNED

    # count leaves
    for node in tree.find_clades(order = 'postorder'):
        node.leaf_count = 1 if node.is_terminal() else np.sum([c.leaf_count for c in node])

    for node in tree.get_nonterminals():
        for c in node:
            c.up=node
    tree.root.up = None
    tree.root.sequences = {'nuc':{}}
    tree.root.sequences.update({gene:{} for gene in all_muts.get(tree.root.name, {}).get('aa_muts', {})})

    # attach sequences to all nodes
    for node in tree.find_clades(order='preorder'):
        if node.up:
            node.sequences = {gene:muts.copy() for gene, muts in node.up.sequences.items()}
        for mut in all_muts.get(node.name, {}).get('muts', []):
            a, pos, d = mut[0], int(mut[1:-1])-1, mut[-1]
            node.sequences['nuc'][pos] = d
        if 'aa_muts' in all_muts.get(node.name, {}):
            for gene in all_muts[node.name]['aa_muts']:
                for mut in all_muts[node.name]['aa_muts'][gene]:
                    a, pos, d = mut[0], int(mut[1:-1])-1, mut[-1]

                    if gene not in node.sequences:
                        node.sequences[gene]={}
                    node.sequences[gene][pos] = d


    # second pass to assign basal nodes within each clade to the clade_labels dict
    # if multiple nodes match, assign annotation to largest
    # otherwise occasional unwanted cousin nodes get assigned the annotation
    for clade_name, clade_alleles in clade_designations.items():
        node_counts = []
        for node in tree.find_clades(order = 'preorder'):
            if is_node_in_clade(clade_alleles, node, ref):
                node_counts.append(node)
        sorted_nodes = sorted(node_counts, key=lambda x: x.leaf_count, reverse=True)
        if len(sorted_nodes) > 0:
            target_node = sorted_nodes[0]
            clade_membership[target_node.name] = clade_name
            clade_labels[target_node.name] = clade_name

    # third pass to propagate clade_membership to descendant nodes
    # (until we encounter a node with its own clade_membership)
    for node in tree.find_clades(order = 'preorder'):
        for child in node:
            if child.name not in clade_labels:
                clade_membership[child.name] = clade_membership[node.name]

    return (clade_membership, clade_labels)

def warn_if_clades_not_found(membership, clade_designations):
    clades = set(clade_designations.keys())
    found = set([clade for clade in membership.values() if clade!=UNASSIGNED])
    if not(len(found)):
        print(f"WARNING in augur.clades: no clades found in tree!")
        return
    for clade in clades-found:
        # warn loudly - one line per unfound clade
        print(f"WARNING in augur.clades: clade '{clade}' not found in tree!")


def get_reference_sequence_from_root_node(all_muts, root_name):
    """
    Extracts the (nuc) sequence from the root node, if set, as well as
    the (aa) sequences. Returns a dictionary of {geneName: rootSequence}
    where rootSequence is a list and geneName may be 'nuc'.
    """
    ref = {}

    # the presence of a single mutation will imply that the corresponding reference
    # sequence should be present, and we will warn if it is not
    nt_present = False
    genes_present = set([])
    missing = []
    for d in all_muts.values():
        if "muts" in d:
            nt_present = True
        genes_present.update(d.get('aa_muts', {}).keys())

    if nt_present:
        try:
            ref['nuc'] = list(all_muts.get(root_name, {})["sequence"])
        except KeyError:
            missing.append("nuc")

    for gene in genes_present:
        try:
            ref[gene] = list(all_muts.get(root_name, {}).get("aa_sequences", {})[gene])
        except KeyError:
            missing.append(gene)

    if missing:
        print(f"WARNING in augur.clades: sequences at the root node have not been specified for {{{', '.join(missing)}}}, \
even though mutations were observed. Clades which are annotated using bases/codons present at the root \
of the tree may not be correctly inferred.")

    return ref

def parse_nodes(tree_file, node_data_files, validation_mode):
    tree = Phylo.read(tree_file, 'newick')
    # don't supply tree to read_node_data as we don't want to require that every node is present in the node_data JSONs
    node_data = read_node_data(node_data_files, validation_mode=validation_mode)
    # node_data files can be parsed without 'nodes' (if they have 'branches')
    if "nodes" not in node_data or len(node_data['nodes'].keys())==0:
        raise AugurError(f"No nodes found in the supplied node data files. Please check {', '.join(node_data_files)}")
    json_nodes = set(node_data["nodes"].keys())
    tree_nodes = set([clade.name for clade in tree.find_clades()])
    if not json_nodes.issubset(tree_nodes):
        raise AugurError(f"The following nodes in the node_data files ({', '.join(node_data_files)}) are not found in the tree ({tree_file}): {', '.join(json_nodes - tree_nodes)}")
    ensure_no_multiple_mutations(node_data['nodes'])
    return (tree, node_data['nodes'])

def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("clades", help=__doc__)
    parser.add_argument('--tree', required=True, help="prebuilt Newick -- no tree will be built if provided")
    parser.add_argument('--mutations', required=True, metavar="NODE_DATA_JSON", nargs='+', action=ExtendOverwriteDefault, help='JSON(s) containing ancestral and tip nucleotide and/or amino-acid mutations ')
    parser.add_argument('--reference', nargs='+', action=ExtendOverwriteDefault, help=SUPPRESS)
    parser.add_argument('--clades', required=True, metavar="TSV", type=str, help='TSV file containing clade definitions by amino-acid')
    parser.add_argument('--output-node-data', type=str,  metavar="NODE_DATA_JSON", help='name of JSON file to save clade assignments to')
    parser.add_argument('--membership-name', type=str, default="clade_membership", help='Key to store clade membership under; use "None" to not export this')
    parser.add_argument('--label-name', type=str, default="clade", help='Key to store clade labels under; use "None" to not export this')
    add_validation_arguments(parser)
    return parser


def run(args):
    (tree, all_muts) = parse_nodes(args.tree, args.mutations, args.validation_mode)

    if args.reference:
        # PLACE HOLDER FOR vcf WORKFLOW.
        # Works without a reference for now but can be added if clade defs contain positions
        # that are monomorphic across reference and sequence sample.
        print(f"WARNING in augur.clades: You have provided a --reference file(s) ({args.reference}) however this is functionality is not yet supported.")
        ref = None
    else:
        # extract reference sequences from the root node entry in the mutation json
        # if this doesn't exist, it will complain but not error.
        ref = get_reference_sequence_from_root_node(all_muts, tree.root.name)

    clade_designations = read_in_clade_definitions(args.clades)
    membership, labels = assign_clades(clade_designations, all_muts, tree, ref)
    warn_if_clades_not_found(membership, clade_designations)

    membership_key= args.membership_name if args.membership_name.upper() != "NONE" else None
    label_key= args.label_name if args.label_name.upper() != "NONE" else None

    node_data_json = {}
    if membership_key:
        node_data_json['nodes'] = {node: {membership_key: clade} for node,clade in membership.items()}
        print(f"Clade membership stored on nodes → <node_name> → {membership_key}", file=sys.stdout)
    if label_key:
        node_data_json['branches'] = {node: {'labels': {label_key: label}} for node,label in labels.items()}
        print(f"Clade labels stored on branches → <node_name> → labels → {label_key}", file=sys.stdout)

    out_name = get_json_name(args)
    write_json(node_data_json, out_name)
    print(f"Clades written to {out_name}", file=sys.stdout)
