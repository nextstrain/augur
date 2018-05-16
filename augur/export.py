import numpy as np
from Bio import Phylo
from collections import defaultdict
from .utils import read_metadata, read_node_data, write_json, read_config

def tree_to_json(node, extra_attr = []):
    '''
    converts the Biopython tree structure to a dictionary that can
    be written to file as a json. This is called recursively.
    input
        node -- node for which top level dict is produced.
        extra_attr -- attributes to export in addition to strain name and numdate
    '''
    tree_json = {}
    str_attr = ['strain']
    num_attr = ['numdate']
    if hasattr(node, 'name'):
        tree_json['strain'] = node.name

    # transfer attributes, round numerical ones
    for prop in str_attr:
        if hasattr(node, prop):
            tree_json[prop] = node.__getattribute__(prop)
    for prop in num_attr:
        if hasattr(node, prop):
            try:
                tree_json[prop] = round(node.__getattribute__(prop),5)
            except:
                print("cannot round:", node.__getattribute__(prop), "assigned as is")
                tree_json[prop] = node.__getattribute__(prop)

    # loop over extra attributes. Extra-attributed be tuples where
    # element 0 is the name of the attribute, 1 is a callable to
    # convert the data to the desired format
    for prop in extra_attr:
        if len(prop)==2 and callable(prop[1]):
            if hasattr(node, prop[0]):
                tree_json[prop] = prop[1](node.__getattribute__(prop[0]))
        else:
            if hasattr(node, prop):
                tree_json[prop] = node.__getattribute__(prop)

    # call on children
    if node.clades:
        tree_json["children"] = []
        for ch in node.clades:
            tree_json["children"].append(tree_to_json(ch, extra_attr))

    return tree_json

# put the data saved in node_data json back onto the tree
# we really should just grab these data during export.
def attach_tree_meta_data(T, node_meta):
    def process_mutations(muts):
        realMut = [a+str(pos+1)+d for (a, pos ,d) in muts if a!='-' and d!='-']
        #This exclude gaps from displaying in Auspice. They're ignored, in
        #treebuilding anyway, and clutter up display/prevent from seeing real ones
        if len(realMut)==0:
            realMut = [""]
        return realMut

    for n in T.find_clades(order='preorder'):
        if n.name not in node_meta:
            print("ERROR: keys in tree and node meta data don't match. Node %s is missing"%n.name)
            continue

        n.attr={}
        n.aa_muts={}
        for field, val in node_meta[n.name].items():
            if field=='sequence':
                continue
            if field=='mutations':
                muts = process_mutations(val)
                if muts[0]: #must test if string empty, not array! causes NaN error in Auspice
                    n.__setattr__('muts', muts)
            elif field=='aa_muts':
                n.aa_muts = {}
                for prot, tmp_muts in val.items():
                    muts = process_mutations(tmp_muts)
                    if muts[0]: #must test if string is empty, not array!
                        n.aa_muts[prot] = muts
            elif field in ['branch_length', 'mutation_length', 'clock_length',
                           'clade', 'num_date', 'numdate']:
                n.__setattr__(field, val)
                n.attr[field] = val
            else:
                n.attr[field] = val
            if field=='numdate':
                n.__setattr__("num_date", val)
                n.attr["num_date"] = val

    # calculate divergence by summing branch length
    T.root.attr['div']=0
    for n in T.get_nonterminals(order='preorder'):
        for c in n:
            bl =  c.mutation_length if hasattr(c, "mutation_length") else c.branch_length
            c.attr["div"] = n.attr["div"] + bl

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
        n.xvalue = n.attr['div']


def summarise_publications(metadata):
    info = defaultdict(lambda: {"n": 0, "title": "?"})
    mapping = {}
    for n, d in metadata.items():
        if "authors" not in d:
            mapping[n] = None
            print("Error - {} had no authors".format(n))
            continue

        authors = d["authors"]
        mapping[n] = authors
        info[authors]["n"] += 1
        for attr in ["title", "journal", "paper_url"]:
            if attr in d:
                info[authors][attr] = d[attr]

    return (info, mapping)

def make_control_json(T, controls):
    controls_json = {}
    for super_cat, fields in controls.items():
        cat_count = {}
        for n in T.get_terminals():
            tmp = cat_count
            for field in fields:
                tmp["name"] = field
                if field in n.attr:
                    cat = n.attr[field]
                else:
                    cat='unknown'
                if cat in tmp:
                    tmp[cat]['count']+=1
                else:
                    tmp[cat] = {'count':1, 'subcats':{}}
                tmp = tmp[cat]['subcats']
        controls_json[super_cat] = cat_count
    return controls_json

def read_color_maps(fname):
    cm = defaultdict(list)
    try:
        with open(fname) as fh:
            for line in fh:
                # line: trait   trait_value     hex_code
                if line.startswith('#'): continue
                fields = line.strip().split()
                if len(fields)!=3: continue
                cm[fields[0]].append((fields[1], fields[2]))
    except IOError:
        print("WARNING: Couldn't open color definitions file {}.".format(fname))

    return cm


def export_metadata_json(T, metadata, tree_meta, config, color_map_file, fname, indent=0):
    meta_json = {}
    terminals = [n.name for n in T.get_terminals()]
    meta_json["virus_count"] = len(terminals)
    meta_subset = {k:v for k,v in metadata.items() if k in terminals}
    color_maps = read_color_maps(color_map_file)

    (author_info, seq_to_author) = summarise_publications(meta_subset)
    meta_json["author_info"] = author_info
    meta_json["seq_author_map"] = seq_to_author

    # join up config color options with those in the input JSONs.
    col_opts = config["color_options"]
    for trait in col_opts:
        if trait in color_maps:
            col_opts[trait]["color_map"] = color_maps[trait]

    meta_json["annotations"] = tree_meta['annotation']

    meta_json.update(config)
    if len(config["controls"]):
        meta_json["controls"] = make_control_json(T, config["controls"])

    write_json(meta_json, fname)


def run(args):
    # load data, process, and write out
    T = Phylo.read(args.tree, 'newick')
    seq_meta, meta_columns = read_metadata(args.metadata)
    tree_meta = read_node_data(args.node_data, traits=args.traits, aa_muts=args.aa_muts)
    attach_tree_meta_data(T, tree_meta["nodes"])
    tree_layout(T)
    fields_to_export = list(list(tree_meta['nodes'].values())[0].keys())\
                       +["clade","tvalue","yvalue", "xvalue", "attr", "muts", "aa_muts"]
    tjson = tree_to_json(T.root, extra_attr=fields_to_export)
    write_json(tjson, args.output)

    export_metadata_json(T, seq_meta, tree_meta, read_config(args.config),
                         args.color_defs, args.meta_output)
