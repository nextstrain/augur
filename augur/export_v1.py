"""
Augur Version 1-specific helpers for exporting JSON files suitable for
visualization with auspice.
"""

import os, sys
import re
import time
import numpy as np
from Bio import Phylo
from argparse import SUPPRESS
from collections import defaultdict
from .utils import read_metadata, read_node_data, write_json, read_config, read_lat_longs, read_colors

def convert_tree_to_json_structure(node, metadata, div=0, strains=None):
    """
    converts the Biopython tree structure to a dictionary that can
    be written to file as a json. This is called recursively.
    Creates the strain property & divergence on each node

    input
        node -- node for which top level dict is produced.
        div  -- cumulative divergence (root = 0)

    returns
        tree in JSON structure
        list of strains
    """
    node_struct = {
        'attr': {"div": div},
        'strain': node.name,
        'clade': node.clade
    }
    for attr in ['branch_length', 'tvalue', 'yvalue', 'xvalue']:
        try:
            node_struct[attr] = node.__getattribute__(attr)
        except AttributeError:
            pass


    if strains is None:
        strains = [node_struct["strain"]]
    else:
        strains.append(node_struct["strain"])

    if node.clades:
        node_struct["children"] = []
        for child in node.clades:
            if 'mutation_length' in metadata[child.name]:
                cdiv = div + metadata[child.name]['mutation_length']
            elif 'branch_length' in metadata[child.name]:
                cdiv = div + metadata[child.name]['branch_length']
            node_struct["children"].append(convert_tree_to_json_structure(child, metadata, div=cdiv, strains=strains)[0])

    return (node_struct, strains)



def recursively_decorate_tree_json_v1_schema(node, node_metadata, decorations):
    """
    This function is deprecated and is used to produce the v1-compatable JSON format

    For given decorations, add information from node_metadata to
    each node in the tree.
    * decorations must have property "key" which is the key used to insert
    into the node and the default key used to access node_metadata
    * if decorations has property "lookup_key", this is used to access
    node meta_data instead
    * if decorations has property "is_attr" (and it's value is True)
    then the result is inserted into node["attr"]

    returns Null
    """
    try:
        metadata = node_metadata[node["strain"]]
        metadata["strain"] = node["strain"]
    except KeyError:
        raise Exception("ERROR: node %s is not found in the node metadata."%n.name)

    for data in decorations:
        val = None
        insert_key = data["key"]
        try:
            if "lookup_key" in data:
                val = metadata[data["lookup_key"]]
            else:
                val = metadata[insert_key]
        except KeyError:
            pass

        if val is not None:
            if "is_attr" in data and data["is_attr"]:
                node["attr"][insert_key] = val
            else:
                if insert_key == 'aa_muts':
                    val = {k:v for k,v in val.items() if len(v) }
                node[insert_key] = val

    if "children" in node:
        for child in node["children"]:
            recursively_decorate_tree_json_v1_schema(child, node_metadata, decorations)


def tree_layout(T):
    """
    calculate tree layout.
    This function is deprecated, and only used for the v1 JSON format
    """
    yval=T.count_terminals()
    clade = 0;
    for n in T.find_clades(order='postorder'):
        n.clade=clade; clade+=1;
        if n.is_terminal():
            n.yvalue=yval
            yval-=1
        else:
            child_yvalues = [c.yvalue for c in n]
            n.yvalue=0.5*(np.min(child_yvalues)+np.max(child_yvalues))


def process_colorings(jsn, color_mapping, nodes=None, node_metadata=None):
    if "color_options" not in jsn:
        print("WARNING: no color options were defined")
        return
    data = jsn["color_options"]

    for trait, options in data.items():
        if "legendTitle" not in options: options["legendTitle"] = trait
        if "menuItem" not in options: options["menuItem"] = trait
        if "key" not in options: options["key"] = trait

        if nodes:
            values_in_tree = {node[trait] for node in nodes.values() if trait in node}
        else:
            values_in_tree = {data["traits"][trait]["value"] for name, data in node_metadata.items() if trait in data['traits']}

        if trait.lower() in color_mapping:
            # remember that the color maps (from the TSV) are in lower case, but this is not how they should be exported
            case_map = {str(val).lower(): val for val in values_in_tree}
            options["color_map"] = [(case_map[m[0]], m[1]) for m in color_mapping[trait.lower()] if m[0] in case_map]

    return data


def process_geographic_info(jsn, lat_long_mapping, node_metadata=None, nodes=None):
    if "geo" not in jsn:
        return {}
    geo = defaultdict(dict)

    traits = jsn["geo"]

    for trait in traits:
        demes_in_tree = {node[trait] for node in nodes.values() if trait in node}

        for deme in demes_in_tree:
            try:
                geo[trait][deme] = lat_long_mapping[(trait.lower(),deme.lower())]
            except KeyError:
                print("Error. {}->{} did not have an associated lat/long value (matching performed in lower case)".format(trait, deme))
    return geo


def process_annotations(node_data):
    # treetime adds "annotations" to node_data
    if "annotations" not in node_data: # if haven't run tree through treetime
        return None
    return node_data["annotations"]

def process_panels(user_panels, meta_json):
    try:
        panels = meta_json["panels"]
    except KeyError:
        panels = ["tree", "map", "entropy"]

    if user_panels is not None and len(user_panels) != 0:
        panels = user_panels

    if "geo" in meta_json:
        geoTraits = meta_json["geo"].keys()
    else:
        geoTraits = []

    if "annotations" in meta_json:
        annotations = meta_json["annotations"].keys()
    else:
        annotations = []

    if "entropy" in panels and len(annotations) == 0:
        panels.remove("entropy")
    if "map" in panels and len(geoTraits) == 0:
        panels.remove("map")

    return panels


def construct_author_info_v1(metadata, tree, nodes):
    """
    author info maps the "authors" property present on tree nodes
    to further information about the paper etc
    """

    authorsInTree = set()
    for node in tree.find_clades(order='postorder'):
        if node.is_terminal and node.name in nodes and "authors" in nodes[node.name]:
            authorsInTree.add(nodes[node.name]["authors"])

    author_info = defaultdict(lambda: {"n": 0})
    for strain, data in metadata.items():
        if "authors" not in data:
            print("Error - {} had no authors".format(strain))
            continue
        if data["authors"] not in authorsInTree:
            continue
        authors = data["authors"]
        author_info[authors]["n"] += 1
        # add in extra attributes if they're present in the meta TSV (for this strain...)
        for attr in ["title", "journal", "paper_url"]:
            if attr in data:
                if attr in author_info[authors] and data[attr].strip() != author_info[authors][attr].strip():
                    print("Error - {} had contradictory {}(s): {} vs {}".format(authors, attr, data[attr], author_info[authors][attr]))
                author_info[authors][attr] = data[attr].strip()

    return author_info


def add_tsv_metadata_to_nodes(nodes, meta_tsv, meta_json, extra_fields=['authors', 'url', 'accession']):
    """
    Only used for v1 schema compatability

    Add the relevent fields from meta_tsv to the nodes
    (both are dictionaries keyed off of strain names)
    * the relevent fields are found by scanning the meta json
    together with the extra_fields param
    """
    fields = [x for x in meta_json["color_options"].keys() if x != "gt"] + extra_fields
    if "geo" in meta_json:
        fields += meta_json["geo"]

    for strain, node in nodes.items():
        if strain not in meta_tsv:
            continue
        for field in fields:
            if field not in node and field in meta_tsv[strain] and meta_tsv[strain][field]:
                node[field] = meta_tsv[strain][field]


def get_root_sequence(root_node, ref=None, translations=None):
    '''
    create a json structure that contains the sequence of the root, both as
    nucleotide and as translations. This allows look-up of the sequence for
    all states, including those that are not variable.

    Parameters
    ----------
    root_node : dict
    	data associated with the node
    ref : str, optional
        filename of the root sequence
    translations : str, optional
        file name of translations

    Returns
    -------
    dict
        dict of nucleotide sequence and translations
    '''
    root_sequence = {}
    if ref and translations:
        from Bio import SeqIO
        refseq = SeqIO.read(ref, 'fasta')
        root_sequence['nuc']=str(refseq.seq)
        for gene in SeqIO.parse(translations, 'fasta'):
            root_sequence[gene.id] = str(gene.seq)
    else:
        root_sequence["nuc"] = root_node["sequence"]
        root_sequence.update(root_node["aa_sequences"])

    return root_sequence


def add_core_args(parser):
    core = parser.add_argument_group("REQUIRED")
    core.add_argument('--tree','-t', required=True, help="tree to perform trait reconstruction on")
    core.add_argument('--metadata', required=True, help="tsv file with sequence meta data")
    core.add_argument('--node-data', required=True, nargs='+', help="JSON files with meta data for each node")
    core.add_argument('--output-tree', help="JSON file name that is passed on to auspice (e.g., zika_tree.json).")
    core.add_argument('--output-meta', help="JSON file name that is passed on to auspice (e.g., zika_meta.json).")
    core.add_argument('--auspice-config', help="file with auspice configuration")
    return core


def add_option_args(parser):
    options = parser.add_argument_group("OPTIONS")
    options.add_argument('--colors', help="file with color definitions")
    options.add_argument('--lat-longs', help="file latitudes and longitudes, overrides built in mappings")
    options.add_argument('--tree-name', default=False, help="Tree name (needed for tangle tree functionality)")
    options.add_argument('--minify-json', action="store_true", help="export JSONs without indentation or line returns")
    options.add_argument('--output-sequence', help="JSON file name that is passed on to auspice (e.g., zika_seq.json).")
    options.add_argument('--reference', required=False, help="reference sequence for export to browser, only vcf")
    options.add_argument('--reference-translations', required=False, help="reference translations for export to browser, only vcf")
    return options


def register_arguments_v1(subparsers):
    # V1 sub-command
    v1 = subparsers.add_parser('v1', help="Export version 1 JSON schema (separate meta and tree JSONs)")
    v1_core = add_core_args(v1)
    v1_options = add_option_args(v1)
    v1.add_argument("--v1", help=SUPPRESS, default=True)


def run_v1(args):
    T = Phylo.read(args.tree, 'newick')
    node_data = read_node_data(args.node_data) # args.node_data is an array of multiple files (or a single file)
    nodes = node_data["nodes"] # this is the per-node metadata produced by various augur modules

    if args.minify_json:
        json_indent = None
    else:
        json_indent = 2

    # export reference sequence data including translations. This is either the
    # inferred sequence of the root, or the reference sequence with respect to
    # which mutations are made on the tree (including possible mutations leading
    # to the root of the tree -- typical case for vcf input data).
    if args.output_sequence:
        if T.root.name in nodes:
            root_sequence = get_root_sequence(nodes[T.root.name], ref=args.reference,
                                                translations=args.reference_translations)
        else:
            root_sequence = {}

        write_json(root_sequence, args.output_sequence)

    meta_json = read_config(args.auspice_config)
    meta_tsv, _ = read_metadata(args.metadata)
    add_tsv_metadata_to_nodes(nodes, meta_tsv, meta_json)

    tree_layout(T)
    tree_json, _ = convert_tree_to_json_structure(T.root, nodes)

    # now the messy bit about what decorations (e.g. "country", "aa_muts") do we want to add to the tree?
    # see recursively_decorate_tree_json to understand the tree_decorations structure
    tree_decorations = [
        {"key": "num_date", "lookup_key": "numdate", "is_attr": True},
        {"key": "muts", "is_attr": False},
        {"key": "aa_muts", "is_attr": False}
    ]
    traits_via_node_metadata = {k for node in nodes.values() for k in node.keys()}
    traits_via_node_metadata -= {'sequence', 'mutation_length', 'branch_length', 'numdate',
                                    'mutations', 'muts', 'aa_muts', 'aa_sequences'}
    for trait in traits_via_node_metadata:
        tree_decorations.append({"key": trait, "is_attr": True})

    recursively_decorate_tree_json_v1_schema(tree_json, nodes, decorations=tree_decorations)
    write_json(tree_json, args.output_tree, indent=json_indent)

    # Export the metadata JSON
    lat_long_mapping = read_lat_longs(args.lat_longs)
    color_mapping = read_colors(args.colors)
    meta_json["updated"] = time.strftime("%d %b %Y")
    meta_json["virus_count"] = len(list(T.get_terminals()))
    meta_json["author_info"] = construct_author_info_v1(meta_tsv, T, nodes)
    meta_json["color_options"] = process_colorings(meta_json, color_mapping, nodes=nodes)
    meta_json["geo"] = process_geographic_info(meta_json, lat_long_mapping, nodes=nodes)
    annotations = process_annotations(node_data)
    if annotations:
        meta_json["annotations"] = annotations
    meta_json["panels"] = process_panels(None, meta_json)

    write_json(meta_json, args.output_meta, indent=json_indent)
    return 0
