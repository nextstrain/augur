import os
import numpy as np
from Bio import Phylo
from collections import defaultdict
from .utils import read_metadata, read_node_data, write_json, read_config, read_lat_longs, read_colors, attach_tree_meta_data

def tree_to_json(node, fields_to_export = [], top_level = [], div=0):
    '''
    converts the Biopython tree structure to a dictionary that can
    be written to file as a json. This is called recursively.
    input
        node -- node for which top level dict is produced.
        fields_to_export -- attributes to export in addition to strain name and numdate
        top_level -- list of strings or tuples of length 2. If tuple, second field is the fn which
                     is called to produce the value. If string, the key->value lookup is used.
    '''
    tree_json = {'attr':{"div":div}, 'branch_length':node.branch_length}
    if hasattr(node, 'name'):
        tree_json['strain'] = node.name

    # transfer attributes, round numerical ones
    for field in fields_to_export+top_level:
        val=None
        if len(field)==2 and callable(field[1]):
            fname = field[0]
            if hasattr(node, fname):
                val = field[1](node.__getattribute__(fname))
        else:
            fname = field
            if hasattr(node, fname):
                val = node.__getattribute__(fname)

        # shadow clade by strain. clade is deprecated and will be removed.
        if field == "clade":
            try: val = node.__getattribute__("strain")
            except: pass

        if field in top_level:
            tree_json[fname] = val
        else:
            tree_json['attr'][fname] = val

    # call on children
    if node.clades:
        tree_json["children"] = []
        for ch in node.clades:
            cdiv = div + (ch.mutation_length if hasattr(ch, "mutation_length") else ch.branch_length)
            tree_json["children"].append(tree_to_json(ch, fields_to_export, top_level, div=cdiv))

    return tree_json

# calculate tree layout. should be obsolete with future auspice versions
def tree_layout(T):
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


def summarise_publications(metadata):
    # Return one format to go into 'author_info' and one to go into 'authors' in 'controls'
    control_authors = defaultdict(lambda: {"count": 0, "subcats": {} })
    author_info = defaultdict(lambda: {"n": 0, "title": "?" })
    for n, d in metadata.items():
        if "authors" not in d:
            print("Error - {} had no authors".format(n))
            continue

        authors = d["authors"]
        author_info[authors]["n"] += 1
        control_authors[authors]["count"] += 1
        for attr in ["title", "journal", "paper_url"]:
            if attr in d:
                author_info[authors][attr] = d[attr]

    return (author_info, control_authors)

def export_metadata_json(T, metadata, tree_meta, config, color_mapping, lat_long_mapping, fname, indent=0):
    meta_json = {}
    import time
    meta_json["updated"] = time.strftime("%d %b %Y")
    terminals = [n.name for n in T.get_terminals()]
    meta_json["virus_count"] = len(terminals)
    meta_subset = {k:v for k,v in metadata.items() if k in terminals}

    (author_info, control_authors) = summarise_publications(meta_subset)
    meta_json["author_info"] = author_info

    if "annotations" not in tree_meta: #if haven't run tree through treetime
        meta_json["annotations"] = {}
        config["panels"] = ["tree","map"]
    if "panels" not in config:
        config["panels"] = ["tree","map","entropy"]

    # join up config color options with those in the input JSONs.
    # TODO: change the schema for these
    col_opts = config["color_options"]
    for trait in col_opts:
        if trait in color_mapping:
            col_opts[trait]["legendTitle"] = trait
            col_opts[trait]["menuItem"] = trait
            col_opts[trait]["key"] = trait
            col_opts[trait]["color_map"] = color_mapping[trait]

    if "annotations" in tree_meta:
        meta_json["annotations"] = tree_meta['annotations']

    meta_json.update(config)

    if "geo" in config:
        geo={}
        for geo_field in config["geo"]:
            geo[geo_field]={}
            for node, attrs in tree_meta["nodes"].items():
                if geo_field in attrs:
                    loc = attrs[geo_field]
                    if loc not in geo[geo_field]:
                        if (geo_field,loc) in lat_long_mapping:
                            geo[geo_field][loc] = lat_long_mapping[(geo_field,loc)]
                        else:
                            print("Lat/long for " + loc + " absent, defaulting to 0,0")
                            geo[geo_field][loc] = {"latitude": 0, "longitude": 0}
        meta_json["geo"] = geo
    write_json(meta_json, fname)


def tree_meta_info(tree_meta, seq_meta, fields=['authors', 'url', 'accession']):
    '''
    Attaches author name to nodes so it's included in tree.json
    Should also perhaps be attaching paper_url, journal, url here?
    '''
    for key in tree_meta['nodes'].keys():
        val = tree_meta['nodes'][key]
        if key in seq_meta:
            for f in fields:
                if f in seq_meta[key]:
                    val[f] = seq_meta[key][f]
    return tree_meta


def run(args):
    # load data, process, and write out
    T = Phylo.read(args.tree, 'newick')
    seq_meta, meta_columns = read_metadata(args.metadata)
    tree_meta = read_node_data(args.node_data) # an array of multiple files (or a single file)
    tree_meta = tree_meta_info(tree_meta, seq_meta) #attach author name to node

    # TODO: remove the following function and combine it with tree_to_json
    attach_tree_meta_data(T, tree_meta["nodes"])

    # TODO: can remove the y-values, they're now calculated in auspice (maybe keep temporarily?)
    tree_layout(T)

    # from the fields in the tree_meta (taken from one or more JSONs), which ones do we want to export?
    node_fields = set()
    for n in tree_meta['nodes'].values():
        node_fields.update(n.keys())
    fields_to_export = [x for x in  node_fields
                        if x not in ['sequence', 'mutations', 'muts', 'aa_muts']]+['num_date']
    # and which fields are top level? (the rest are all in the attr dict)
    top_level = ["clade","tvalue","yvalue", "xvalue", "muts", "aa_muts"]

    tree_json = tree_to_json(T.root, fields_to_export=fields_to_export, top_level=top_level)
    write_json(tree_json, args.output_tree)

    # load defaults and supp file that overrides defaults
    lat_long_mapping = read_lat_longs(args.lat_longs)
    color_mapping = read_colors(args.colors)

    export_metadata_json(T, seq_meta, tree_meta, read_config(args.auspice_config),
                         color_mapping, lat_long_mapping, args.output_meta)
