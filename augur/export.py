import numpy as np
from Bio import Phylo
from .utils import read_metadata, read_nodedata, write_json

def tree_to_json(node, extra_attr = []):
    tree_json = {}
    str_attr = ['strain']
    num_attr = ['numdate']
    if hasattr(node, 'name'):
        tree_json['strain'] = node.name

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

    for prop in extra_attr:
        if len(prop)==2 and callable(prop[1]):
            if hasattr(node, prop[0]):
                tree_json[prop] = prop[1](node.__getattribute__(prop[0]))
        else:
            if hasattr(node, prop):
                tree_json[prop] = node.__getattribute__(prop)

    tree_json['tvalue'] = tree_json['numdate']
    if node.clades:
        tree_json["children"] = []
        for ch in node.clades:
            tree_json["children"].append(tree_to_json(ch, extra_attr))
    return tree_json

def attach_tree_meta_data(T, node_meta):
    def parse_mutations(muts):
        allMut = muts.split(',') if type(muts) in [str] else ""
        realMut=[]
        #This exclude gaps from displaying in Auspice. They're ignored, in
        #treebuilding anyway, and clutter up display/prevent from seeing real ones
        for m in allMut:
            if '-' not in m:
                realMut.append(m)
        if len(realMut)==0:
            realMut = [""]
        return realMut
        #return muts.split(',') if type(muts) in [str, unicode] else ""

    for n in T.find_clades(order='preorder'):
        n.attr={}
        n.aa_muts={}
        for field, val in node_meta[n.name].items():
            if 'mutations' in field:
                if field=='mutations':
                    muts = parse_mutations(val)
                    if muts[0]: #must test if string empty, not array! causes NaN error in Auspice
                        n.__setattr__('muts', muts)
                else:
                    prot = '_'.join(field.split('_')[:-1])
                    muts = parse_mutations(val)
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

    T.root.attr['div']=0
    for n in T.get_nonterminals(order='preorder'):
        for c in n:
            bl =  c.mutation_length if hasattr(c, "mutation_length") else c.branch_length
            c.attr["div"] = n.attr["div"] + bl


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


def run(args):
    T = Phylo.read(args.tree, 'newick')
    tree_meta = read_nodedata(args.nodedata, traits=args.traits)['nodes']
    attach_tree_meta_data(T, tree_meta)
    tree_layout(T)
    fields_to_export = list(list(tree_meta.values())[0].keys())\
                       +["clade","tvalue","yvalue", "xvalue", "attr", "muts", "aa_muts"]
    tjson = tree_to_json(T.root, extra_attr=fields_to_export)
    write_json(tjson, args.output)

