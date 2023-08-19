from __future__ import print_function
from Bio import Phylo
import imp
from StringIO import StringIO
import json

def modified_tree_to_json(node, extra_attr = []):
    tree_json = {}
    str_attr = ['strain','attr']
    num_attr = ['xvalue', 'yvalue', 'tvalue', 'num_date', 'clade']
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

    if node.clades:
        tree_json["children"] = []
        for ch in node.clades:
            tree_json["children"].append(modified_tree_to_json(ch, extra_attr))
    return tree_json

def mock_meta_json(n):
    meta = {}
    meta["updated"] = "today"
    meta["virus_count"] = n
    meta["title"] = "BEAST proof of principle (MERS-CoV)"

    meta["color_options"] = {
      "host": {
        "menuItem": "host",
        "type": "discrete",
        "legendTitle": "host",
        "key": "host"
      },
      "type.prob": {
        "menuItem": "host.prob",
        "type": "continuous",
        "legendTitle": "host.prob",
        "key": "type.prob"
      },
        "posterior": {
          "menuItem": "posterior",
          "type": "continuous",
          "legendTitle": "posterior",
          "key": "posterior"
        },
    }
    meta["filters"] = ["type"]
    meta["commit"] = "unknown"
    meta["panels"] = ["tree"]
    meta["geo"] = {}
    return meta;


if __name__=="__main__":
    print("This is a proof of principle script - it should not be relied upon for analysis.")
    try:
        bt = imp.load_source('baltic', '/Users/james/blab/baltic/baltic.py')
    except IOError:
        bt = imp.load_source('baltic', '/Users/evogytis/Documents/BLAB_baltic/baltic.py')

    # set_trace()
    bt_tree=bt.loadNexus('./scratch/MERS_274_sCoal.combinedTyped.mcc.tree',tip_regex='\|([0-9\-]+)$') ## loads a BEAST nexus file
    mostRecentTip=max([k.absoluteTime for k in bt_tree.Objects])

    simplified_tree=bt_tree.toString() ## returns a newick string without beast annotations
    xml_tree=Phylo.read(StringIO(simplified_tree),'newick') ## load simplified newick as a biopython XML tree

    xml_nodes=xml_tree.get_nonterminals(order='preorder') ## fetch nodes from XML tree
    bt_nodes=[k for k in bt_tree.traverse_tree(include_all=True) if k.branchType=='node'] ## fetch nodes from baltic tree

    xml_leaves=xml_tree.get_terminals(order='preorder') ## fetch leaves from XML tree
    bt_leaves=[k for k in bt_tree.traverse_tree()] ## fetch leaves from baltic tree


for pairs in [[xml_nodes,bt_nodes],[xml_leaves,bt_leaves]]:
    for x,b in zip(*pairs): ## iterate over pairs of branches in both trees


        setattr(x,'tvalue',b.height) ## height in years
        setattr(x,'xvalue',b.height) ## 0.0 in the future?
        setattr(x,'yvalue',b.y) ## y value after drawing
        setattr(x,'clade',b.index) ## clade

        attrs={}

        if b.absoluteTime!=None:
            attrs['num_date']=b.absoluteTime
        else:
            attrs['num_date']=bt_tree.Objects[0].absoluteTime

        if b.traits.has_key('height_95%_HPD'):
            attrs['num_date_confidence']=[mostRecentTip-t for t in b.traits['height_95%_HPD']]

        if b.branchType=='node':
            setattr(x,'name',b.index)

        if b.traits.has_key('type'):
            attrs['host'] = b.traits['type']

        for trait in b.traits: ## iterate over MCC traits, transfer into attrs dict
            if 'set.prob' in trait: ## all posterior sets converted to confidences
                tr=trait.split('.')[0]
                confidence={t:p for t,p in zip(b.traits['%s.set'%(tr)],b.traits['%s.set.prob'%(tr)])}
                attrs['%s_confidence'%(tr)]=confidence

            attrs[trait]=b.traits[trait]

        setattr(x,'attr',attrs)

    tree_json = modified_tree_to_json(xml_tree.clade)
    tree_json = tree_json["children"][0]
    json.dump(tree_json, open("scratch/mers_tree.json", 'w'), indent=2)
    meta_json = mock_meta_json(len(bt_leaves))
    json.dump(meta_json, open("scratch/mers_meta.json", 'w'), indent=2)
    json.dump({}, open("scratch/mers_sequences.json", 'w'), indent=2)

    # set_trace()

    # attrs: tree_json['children'][0]['attr'].keys()
