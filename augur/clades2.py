"""
Infer clades at internal nodes from metadata
"""

from collections import defaultdict
import sys

from Bio import Phylo
import numpy as np
from treetime.wrappers import reconstruct_discrete_traits

from .errors import AugurError
from .io.metadata import DEFAULT_DELIMITERS, DEFAULT_ID_COLUMNS, InvalidDelimiter, read_metadata
from .utils import write_json

def mugration_inference(tree=None, seq_meta=None, field=None, missing='?'):
    """
    Infer likely ancestral states of a discrete character assuming a time reversible model.

    Parameters
    ----------
    tree : str
        name of tree file
    seq_meta : pandas.DataFrame
        meta data associated with sequences
    field : str, optional
        meta data field to use
    missing : str, optional
        character that is to be interpreted as missing data, default='?'

    Returns
    -------
    T : Bio.Phylo.BaseTree.Tree
        Biophyton tree
    """

    T = Phylo.read(tree, 'newick')
    traits = {}
    nodes = {n.name:n for n in T.get_terminals()}
    for name, meta in seq_meta.iterrows():
        if field in meta and name in nodes and meta[field] != missing:
            traits[name] = meta[field]
    unique_states = list(set(traits.values()))

    if len(unique_states)==0:
        print("WARNING: no states found for clade reconstruction. Using `Unassigned`.", file=sys.stderr)
        for node in T.find_clades():
            node.__setattr__(field, "Unassigned")
        return T
    elif len(unique_states)==1:
        print("WARNING: only one state found for clade reconstruction:", unique_states, file=sys.stderr)
        for node in T.find_clades():
            node.__setattr__(field, unique_states[0])
        return T
    elif len(unique_states)<300:
        tt, letter_to_state, reverse_alphabet = reconstruct_discrete_traits(T, traits, missing_data=missing)
    else:
        raise AugurError("ERROR: 300 or more distinct discrete states found. TreeTime is currently not set up to handle that many states.")

    if tt is None:
        raise AugurError("ERROR in discrete state reconstruction in TreeTime. Please look for errors above.")

    # attach inferred states as e.g. node.region = 'africa'
    for node in tt.tree.find_clades():
        node.__setattr__(field, letter_to_state[node.cseq[0]])

    return tt.tree


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("clades2", help=__doc__)
    parser.add_argument('--tree', '-t', required=True, help="tree to assign clades on")
    parser.add_argument('--metadata', required=True, metavar="FILE", help="table with metadata")
    parser.add_argument('--metadata-delimiters', default=DEFAULT_DELIMITERS, nargs="+",
                        help="delimiters to accept when reading a metadata file. Only one delimiter will be inferred.")
    parser.add_argument('--metadata-id-columns', default=DEFAULT_ID_COLUMNS, nargs="+",
                        help="names of possible metadata columns containing identifier information, ordered by priority. Only one ID column will be inferred.")
    parser.add_argument('--clade-column', required=True, help='column name of metadata field to use for clade inference')
    parser.add_argument('--output-field-name', type=str, default='clade_membership', help='name of field to save clade inferences to')
    parser.add_argument('--output-node-data', required=True, type=str, help='name of JSON file to save clade inferences to')
    parser.epilog = "Note that missing data must be represented by a `?` character. Missing data will currently be inferred."
    return parser


def run(args):
    """run clade inference

    Parameters
    ----------
    args : argparse.Namespace
        command line arguments are parsed by argparse
    """
    tree_fname = args.tree
    try:
        traits = read_metadata(
            args.metadata,
            delimiters=args.metadata_delimiters,
            id_columns=args.metadata_id_columns)
    except InvalidDelimiter:
        raise AugurError(
                f"Could not determine the delimiter of {args.metadata!r}. "
                f"Valid delimiters are: {args.metadata_delimiters!r}. "
                "This can be changed with --metadata-delimiters."
            )

    T = Phylo.read(tree_fname, 'newick')
    missing_internal_node_names = [n.name is None for n in T.get_nonterminals()]
    if np.all(missing_internal_node_names):
        print("\n*** WARNING: Tree has no internal node names!", file=sys.stderr)
        print("*** Without internal node names, ancestral traits can't be linked up to the correct node later.", file=sys.stderr)
        print("*** If you want to use 'augur export' later, re-run this command with the output of 'augur refine'.", file=sys.stderr)
        print("*** If you haven't run 'augur refine', you can add node names to your tree by running:", file=sys.stderr)
        print("*** augur refine --tree %s --output-tree <filename>.nwk" % (tree_fname), file=sys.stderr)
        print("*** And use <filename>.nwk as the tree when running 'ancestral', 'translate', and 'traits'", file=sys.stderr)


    mugration_states = defaultdict(dict)

    from treetime import version as treetime_version
    print(f"augur clades2 is using TreeTime version {treetime_version}")

    T= mugration_inference(tree=tree_fname, seq_meta=traits, field=args.clade_column)

    if T is None: # something went wrong
        raise AugurError("TreeTime failed to infer ancestral states. Something went wrong.")

    for node in T.find_clades():
        mugration_states[node.name][args.output_field_name] = getattr(node, args.clade_column)


    write_json({"nodes":mugration_states},args.output_node_data)

    print("\nInferred ancestral states of clade character using TreeTime:"
          "\n\tSagulenko et al. TreeTime: Maximum-likelihood phylodynamic analysis"
          "\n\tVirus Evolution, vol 4, https://academic.oup.com/ve/article/4/1/vex042/4794731\n", file=sys.stdout)

    print("results written to", args.output_node_data, file=sys.stdout)
