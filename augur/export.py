import numpy as np
from Bio import Phylo


def tree_to_json(node, extra_attr = []):
    tree_json = {}
    str_attr = ['strain']
    num_attr = ['num_date']
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

    tree_json['tvalue'] = tree_json['num_date']
    if node.clades:
        tree_json["children"] = []
        for ch in node.clades:
            tree_json["children"].append(tree_to_json(ch, extra_attr))
    return tree_json


def attach_tree_meta_data(T, node_meta):
    def parse_mutations(muts):
        allMut = muts.split(',') if type(muts) in [str, unicode] else ""
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
                           'clade', 'num_date']:
                n.__setattr__(field, val)
                n.attr[field] = val
            else:
                n.attr[field] = val

    T.root.attr['div']=0
    for n in T.get_nonterminals(order='preorder'):
        for c in n:
            bl =  c.mutation_length if hasattr(c, "mutation_length") else c.branch_length
            c.attr["div"] = n.attr["div"] + bl


def export_metadata_json(T, path, prefix, reference, isvcf=False, indent=1):
    print("Writing out metaprocess")
    mjson = {}

    mjson["virus_count"] = T.count_terminals()
    from datetime import date
    mjson["updated"] = date.today().strftime('%Y-%m-%d')
    mjson["author_info"] = {
        "?": {
           "paper_url": "?",
           "journal": "?",
           "title": "?",
           "n": 1
        }}
    mjson["seq_author_map"] = {}

    from collections import defaultdict
    cmaps = defaultdict(list)
    with open(color_maps(path), 'r') as cfile:
        for line in cfile:
            try:
                trait, name, color = line.strip().split('\t')
            except:
                continue
            cmaps[trait].append((name, color))

    #if drug-resistance colours have been auto-generated, get these too
    import os.path
    if os.path.isfile(drm_color_maps(path)):
        with open(drm_color_maps(path), 'r') as cfile:
            for line in cfile:
                try:
                    trait, name, color = line.strip().split('\t')
                except:
                    continue
                cmaps[trait].append((name, color))

    mjson["color_options"] = {
      "gt": {
           "menuItem": "genotype",
           "type": "discrete",
           "legendTitle": "Genotype",
           "key": "genotype"
          },
       "num_date": {
           "menuItem": "date",
           "type": "continuous",
           "legendTitle": "Sampling date",
           "key": "num_date"
          }}
    for trait in cmaps:
        mjson["color_options"][trait] = {
        "menuItem":trait,
        "type":"discrete",
        "color_map":cmaps[trait],
        "legendTitle":trait,
        "key":trait
        }

    mjson["panels"] = [
        "tree",
        "map",
        "entropy"
        ]
    mjson["title"] = "NextTB"
    mjson["maintainer"] = "Emma Hodcroft"
    mjson["geo"] = {}
    lat_long_defs = load_lat_long_defs()
    for geo_trait in ['region', "country", 'division']:
        mjson["geo"][geo_trait] = {}
        for n in T.find_clades():
            if geo_trait in n.attr:
                place = n.attr[geo_trait]
                if  (place not in mjson["geo"][geo_trait]
                     and place in lat_long_defs):
                    mjson["geo"][geo_trait][place] = lat_long_defs[place]

    mjson["commit"] = "unknown"
    mjson["filters"] = ["country", "region", "division"]

    genes = load_features(reference)
    anno = {}
    for feat, aln_fname in get_genes_and_alignments(path, tree=False):
        if feat in genes:
            anno[feat] = {"start":int(genes[feat].location.start),
                          "end":int(genes[feat].location.end),
                          "strand":genes[feat].location.strand}

    if isvcf:
        #if vcf, there is no 'gene' called 'nuc' that will be read in
        #above, so manually include it here.
        from filenames import ref_fasta
        from Bio import SeqIO
        refSeq = SeqIO.parse(ref_fasta(path), format='fasta').next()
        anno['nuc'] =   {"start":1,
                         "end":len(refSeq.seq),
                         "strand":1 }

    mjson["annotations"] = anno
    write_json(mjson, meta_json(path,prefix), indent=indent)

def tree_layout(T):
    yval=T.count_terminals()
    for n in T.find_clades(order='postorder'):
        if n.is_terminal():
            n.yvalue=yval
            yval-=1
        else:
            child_yvalues = [c.yvalue for c in n]
            n.yvalue=0.5*(np.min(child_yvalues)+np.max(child_yvalues))
        n.xvalue = n.attr['div']


def run(args):
    T = Phylo.read(args.tree, 'newick')
    seq_meta = read_sequence_meta_data(path)
    tree_meta = read_tree_meta_data(path)
    attach_tree_meta_data(T, tree_meta)
    tree_layout(T)
    fields_to_export = tree_meta.values()[0].keys()+["tvalue","yvalue", "xvalue", "attr", "muts", "aa_muts"]
    tjson = tree_to_json(T.root, extra_attr=fields_to_export)
    write_json(tjson, tree_json(path, args.prefix))

    if not args.vcf:
        #Aupice will soon not need these, and getting them from VCF
        #files is not efficient (processing or memory), so simply
        #do not make if VCF file format.
        export_sequence_json(T, path, args.prefix, indent=1)
        export_diversity(path, args.prefix, args.reference)

    export_metadata_json(T, path, args.prefix, args.reference, args.vcf)

