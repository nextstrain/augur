import os
import time
import numpy as np
from Bio import Phylo
from collections import defaultdict
from .utils import read_metadata, read_node_data, write_json, read_config, read_lat_longs, read_colors

def convert_tree_to_json_structure(node, metadata, div=0):
    """
    converts the Biopython tree structure to a dictionary that can
    be written to file as a json. This is called recursively.
    input
        node -- node for which top level dict is produced.
        div  -- cumulative divergence (root = 0)
    """
    node_struct = {
        'attr': {"div": div},
        'strain': node.name
    }

    # the following are DEPRECATED and to be removed
    for attr in ['branch_length', 'tvalue', 'yvalue', 'xvalue']:
        try:
            node_struct[attr] = node.__getattribute__(attr)
        except AttributeError:
            pass

    if node.clades:
        node_struct["children"] = []
        for child in node.clades:
            cdiv = div + metadata[child.name]['mutation_length']
            node_struct["children"].append(convert_tree_to_json_structure(child, metadata, div=cdiv))

    return node_struct


def recursively_decorate_tree_json(node, node_metadata, decorations):
    """
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
                node[insert_key] = val

    if "children" in node:
        for child in node["children"]:
            recursively_decorate_tree_json(child, node_metadata, decorations)


def tree_layout(T):
    """
    calculate tree layout. should be obsolete with future auspice versions
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


def process_color_options(color_options, color_mapping, nodes):
    for trait, options in color_options.items():
        if "legendTitle" not in options:
            options["legendTitle"] = trait
        if "menuItem" not in options:
            options["menuItem"] = trait
        if "key" not in options:
            options["key"] = trait

        if trait in color_mapping:
            valuesInTree = {node[trait] for node in nodes.values() if trait in node}
            options["color_map"] = [m for m in color_mapping[trait] if m[0] in valuesInTree]

    return color_options

def process_geo_resolutions(meta_json, lat_long_mapping, nodes):
    geo = defaultdict(dict)
    if "geo" not in meta_json:
        return geo
    traits = meta_json["geo"]
    for trait in traits:
        demesInTree = {node[trait] for node in nodes.values() if trait in node}
        for deme in demesInTree:
            try:
                geo[trait][deme] = lat_long_mapping[(trait,deme)]
            except KeyError:
                print("Error. {}->{} did not have an associated lat/long value".format(trait, deme))
    return geo

def process_annotations(node_data):
    # treetime adds "annotations" to node_data
    if "annotations" not in node_data: #if haven't run tree through treetime
        return {}
    return node_data["annotations"]

def process_panels(meta_json):
    try:
        panels = meta_json["panels"]
    except KeyError:
        panels = ["tree", "map", "entropy"]
    if "entropy" in panels and len(meta_json["annotations"].keys()) == 0:
        panels.remove("entropy")
    if "map" in panels and len(meta_json["geo"].keys()) == 0:
        panels.remove("map")
    return panels

def add_tsv_metadata_to_nodes(nodes, meta_tsv, meta_json, extra_fields=['authors', 'url', 'accession']):
    """
    Add the relevent fields from meta_tsv to the nodes
    (both are dictionaries keyed off of strain names)
    * the relevent fields are found by scanning the meta json
    together with the extra_fields param
    """
    fields = [x for x in meta_json["color_options"].keys() if x != "gt"] + meta_json["geo"] + extra_fields

    for strain, node in nodes.items():
        if strain not in meta_tsv:
            continue
        for field in fields:
            if field not in node and field in meta_tsv[strain]:
                node[field] = meta_tsv[strain][field]

def getTerminalKeyValuesFromNodes(tree, nodes, key):
    """
    Find the values for the given key across the tree
    nodes is the per-node metadata dict
    """
    vals = set()
    for node in tree.find_clades(order='postorder'):
        if node.is_terminal and node.name in nodes and key in nodes[node.name]:
            vals.add(nodes[node.name][key])
    return vals

def construct_author_info(metadata, authorsInTree):
    """
    author info maps the "authors" property present on tree nodes
    to further information about the paper etc
    """
    author_info = defaultdict(lambda: {"n": 0})
    for strain, data in metadata.items():
        if "authors" not in data:
            print("Error - {} had no authors".format(n))
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


def run(args):
    T = Phylo.read(args.tree, 'newick')
    meta_tsv, _ = read_metadata(args.metadata)
    node_data = read_node_data(args.node_data) # args.node_data is an array of multiple files (or a single file)
    meta_json = read_config(args.auspice_config)
    nodes = node_data["nodes"] # this is the per-node metadata produced by various augur modules

    # export the tree JSON first
    add_tsv_metadata_to_nodes(nodes, meta_tsv, meta_json)
    tree_layout(T) # TODO: deprecated. Should remove.
    tree_json = convert_tree_to_json_structure(T.root, nodes)

    # now the messy bit about what decorations (e.g. "country", "aa_muts") do we want to add to the tree?
    # see recursively_decorate_tree_json to understand the tree_decorations structure
    tree_decorations = [
        {"key": "clade", "lookup_key": "strain"}, # DEPRECATED. Auspice still refers to "clade" so this must stay for the moment.
        {"key": "num_date", "lookup_key": "numdate", "is_attr": True},
        {"key": "muts", "is_attr": False},
        {"key": "aa_muts", "is_attr": False}
    ]
    traits_via_node_metadata = {k for node in nodes.values() for k in node.keys()}
    traits_via_node_metadata -= {'sequence', 'mutation_length', 'branch_length', 'numdate', 'mutations', 'muts', 'aa_muts', 'clock_length'}
    for trait in traits_via_node_metadata:
        tree_decorations.append({"key": trait, "is_attr": True})

    recursively_decorate_tree_json(tree_json, nodes, decorations=tree_decorations)
    write_json(tree_json, args.output_tree, indent=2)

    # Export the metadata JSON
    lat_long_mapping = read_lat_longs(args.lat_longs)
    color_mapping = read_colors(args.colors)
    meta_json["updated"] = time.strftime("%d %b %Y")
    meta_json["virus_count"] = len(list(T.get_terminals()))
    meta_json["author_info"] = construct_author_info(meta_tsv, getTerminalKeyValuesFromNodes(T, nodes, "authors"))
    meta_json["color_options"] = process_color_options(meta_json["color_options"], color_mapping, nodes)
    meta_json["geo"] = process_geo_resolutions(meta_json, lat_long_mapping, nodes)
    meta_json["annotations"] = process_annotations(node_data)
    meta_json["panels"] = process_panels(meta_json)

    write_json(meta_json, args.output_meta, indent=2)
