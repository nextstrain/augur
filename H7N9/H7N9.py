from __future__ import division, print_function
from collections import defaultdict
import sys
sys.path.append('')  # need to import from base
# sys.path.append('/home/richard/Projects')  # need to import from base
from base.io_util import make_dir, remove_dir, tree_to_json, write_json, myopen
from base.sequences import sequence_set, num_date
from base.tree import tree
from base.process import process
from base.frequencies import alignment_frequencies, tree_frequencies
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import numpy as np
from datetime import datetime
from pdb import set_trace
from pprint import pprint
import matplotlib as mpl
from matplotlib import pyplot as plt
from Bio.Align import MultipleSeqAlignment

"""
This script is modelled on zika.py, but attempts to be more modular
"""

def generate_cmap(data, discrete):
    '''
    data is set (e.g. countries), returns list of list with country -> hex
    '''
    norm = mpl.colors.Normalize(0, len(data) - 1)
    cmap = mpl.cm.get_cmap("viridis")
    if discrete:
        if len(data) <= 10:
            cmap = mpl.cm.get_cmap("Vega10")
        elif len(data) <= 20:
            cmap = mpl.cm.get_cmap("Vega20")

    ret = []
    for idx, val in enumerate(list(data)):
        ret.append([val, mpl.colors.to_hex(cmap(norm(idx)))])
    return ret

def getAllAttrs(key, segments):
    '''
    key = country (e.g.), get set of all countries for the key
    '''
    ret = [x.attributes[key] for x in segments[0].seqs.all_seqs.values()]
    for idx in xrange(1, len(segments)):
        ret += [x.attributes[key] for x in segments[idx].seqs.all_seqs.values()]
    return set(ret)

def save_as_nexus(tree, fname):
    def format_string_attr(node, key):
        return "{}=\"{}\"".format(key, node["attr"][key])

    def stringify_node(node, prev_div, terminal):
        if terminal:
            taxa.append(node["strain"])
        extra = [format_string_attr(node, x) for x in ["country", "host", "division"]]
        return "{}[&{}]:{}".format(len(taxa) if terminal else "", ",".join(extra), float(node["attr"]["div"]) - float(prev_div))

    def tree_walk(node, prev_div):
        if "children" in node:
            subtrees = ",".join([tree_walk(child, node["attr"]["div"]) for child in node["children"]])
            return "({}){}".format(subtrees,stringify_node(node, prev_div, False))
        else:
            return stringify_node(node, prev_div, True)

    def nexus(tree):
        nex = []
        nex += ["#NEXUS", ""]
        nex += ["Begin taxa;", "\tDimensions ntax={};".format(len(taxa)), "\tTaxlabels"]
        nex += ["\t\t"+name for name in taxa]
        nex += ["\t\t;", "End;", ""]
        nex += ["Begin trees;", "\tTranslate"]
        nex += ["\t\t{} {},".format(idx+1, name) for idx, name in enumerate(taxa)]
        nex[-1] = nex[-1][:-1] # remove the final comma
        nex += ["\t\t;", "tree TREE1 = [&R] "+tree+";"]
        nex += ["End;"]
        return nex

    print("Saving to nexus tree {}".format(fname))
    taxa = [] # the id of a name is its (0 based) idx + 1, populated by tree_walk()
    json = tree_to_json(tree.root , ['clade', 'attr']);
    nex = nexus(tree_walk(json, 0))
    with open(fname, 'w') as f:
        f.write("\n".join(nex))

if __name__=="__main__":
    import argparse

    # PARAMETERS:
    params = {
        "HA": {
            "lineage": "flu_H7N9_HA",
            "reference_fname": "H7N9/reference_segments/ha.gb",
            "input_data": "../fauna/data/h7n9_ha",
            "proteins": ['SigPep', 'HA1', 'HA2'],
            "min_seq_length": 1500,
        },
        "NA": {
            "lineage": "flu_H7N9_NA",
            "reference_fname": "H7N9/reference_segments/na.gb",
            "input_data": "../fauna/data/h7n9_na",
            "proteins": ['NA'],
            "min_seq_length": 1200,
        },
        # >A/chicken/Jiangxi/10877/2014|    h7n9|   EPI594171|      2014-02-18|         chicken|    china|      china|      china|      china|          egg|            other_database_import
        #        0                           1         2               3                   4           5           6           7            8          9                   10
        #       'strain',                  'virus', 'accession', 'collection_date',    'host',     'region',   'country', 'division', 'location', 'passage_category', 'submitting_lab']
        "fasta_fields": {
            0:'strain', 2:'accession', 3:'date', 4:'host', 6:'country', 7: 'division'
        },
        "dropped_strains": [
            "A/Chicken/Netherlands/16007311-037041/2016", # not part of epi clade
            "A/Chicken/Netherlands/1600", # not part of epi clade
            "A/duck/Zhejiang/LS02/2014", # not of part of epi clade
            "A/BritishColumbia/1/2015" # travel case. throws off map
        ],
        "viruses_per_month": 500,
        # "viruses_per_month": 10,
        "geo_inference": ['country', 'division', 'host'],
        "confidence": False,
        "earliest_sample": datetime(2013,1,1).date()
        # "earliest_sample": datetime(2016,1,1).date()
    }


    # AUSPICE STUFF
    auspice = {
        "panels": ['tree', 'map', 'entropy'],
        "controls": {'geographic location':['country'], 'authors':['authors']},
        "color_options": {
            "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete", "color_map": []},
            "division":{"key":"division", "legendTitle":"Division", "menuItem":"division", "type":"discrete", "color_map": []},
            "host":{"key":"host", "legendTitle":"Host", "menuItem":"host", "type":"discrete", "color_map": []},
            "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
            "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"}
        }
    }

    HA = process(input_data_path = params["HA"]["input_data"],
       store_data_path = 'store/' + params["HA"]["lineage"] + '_',
       build_data_path = 'build/' + params["HA"]["lineage"] + '_',
       reference=params["HA"]["reference_fname"],
       lat_long_fname='../fauna/source-data/geo_lat_long.tsv',
       proteins=params["HA"]['proteins'],
       method='SLSQP', verbose=0)
    NA = process(input_data_path = params["NA"]["input_data"],
       store_data_path = 'store/' + params["NA"]["lineage"] + '_',
       build_data_path = 'build/' + params["NA"]["lineage"] + '_',
       reference=params["NA"]["reference_fname"],
       lat_long_fname='../fauna/source-data/geo_lat_long.tsv',
       proteins=params["NA"]['proteins'],
       method='SLSQP', verbose=0)

    segments = [HA, NA]
    segmentNames = ["HA", "NA"]

    # STUFF IN COMMON WITH BOTH:
    for idx, segment in enumerate(segments):
        print("SEGMENT: {}".format(segmentNames[idx]))
        segment.load_sequences(fields=params["fasta_fields"])
        # segment.seqs.filter(lambda s: s.attributes["country"] != "?")
        # segment.seqs.filter(lambda s: s.attributes["division"] != "?")
        segment.seqs.filter(lambda s: s.attributes["host"] != "laboratoryderived")
        segment.seqs.filter(lambda s: s.attributes["host"] != "watersample")
        for key in segment.seqs.all_seqs:
            if segment.seqs.all_seqs[key].attributes["division"] == "china":
                    segment.seqs.all_seqs[key].attributes["division"] = "?"

        segment.seqs.filter(lambda s: s.attributes['date']>=params["earliest_sample"])
        segment.seqs.filter(lambda s: s.id not in params["dropped_strains"])
        segment.seqs.filter(lambda s: len(s.seq)>=params[segmentNames[idx]]["min_seq_length"])
        segment.seqs.subsample(category = lambda x:(x.attributes['date'].year, x.attributes['date'].month), threshold=params["viruses_per_month"])

    # HOTSWAP IN THE CMAP DATA NOW THAT WE KNOW WHAT WE'VE GOT
    auspice["color_options"]["country"]["color_map"] = generate_cmap(getAllAttrs("country", segments), True)
    auspice["color_options"]["division"]["color_map"] = generate_cmap(getAllAttrs("division", segments), True)
    auspice["color_options"]["host"]["color_map"] = generate_cmap(getAllAttrs("host", segments), True)

    for idx, segment in enumerate(segments):
        print("SEGMENT: {}".format(segmentNames[idx]))
        # should save un-aligned MFA here
        segment.align(codon_align=False)
        SeqIO.write(segment.seqs.aln, segment.build_data_path + "aligned.mfa", "fasta")
        segment.build_tree(num_distinct_starting_trees=1)
        segment.clock_filter(n_iqd=3, plot=True)
        segment.annotate_tree(Tc=0.02, timetree=True, reroot='best', confidence=params["confidence"])
        for geo_attr in params["geo_inference"]:
            print("running geo inference for {}".format(geo_attr))
            segment.tree.geo_inference(geo_attr)
        # SAVE THE TREE IN NEWICK FORMAT
        # from Bio import Phylo
        # Phylo.write(segment.tree.tree, segment.build_data_path + "timeTree.new", 'newick')
        # SAVE THE TREE IN NEXUS FORMAT
        save_as_nexus(segment.tree.tree, segment.build_data_path + "timeTree.nex")


    # determine ladder rank strains in of every tree
    ladder_ranks = defaultdict(list)
    for seg in segments:
        for leaf in seg.tree.tree.get_terminals():
            ladder_ranks[leaf.name].append(leaf.yvalue)
    for seg in segments:
        for leaf in seg.tree.tree.get_terminals():
            leaf.attr['ladder_ranks'] = np.array(ladder_ranks[leaf.name])
            # print(leaf.attr['ladder_ranks'])
            leaf.nleafs=1
    for seg in segments:
        for node in seg.tree.tree.get_nonterminals(order='postorder'):
            node.nleafs = np.sum([x.nleafs for x in node])
            node.attr['ladder_ranks'] = np.sum([x.nleafs*x.attr['ladder_ranks'] for x in node], axis=0)/node.nleafs
    for seg in segments:
        for leaf in seg.tree.tree.find_clades():
            leaf.attr['ladder_ranks'] = list(leaf.attr['ladder_ranks'])


    for idx, segment in enumerate(segments):
        segment.export(controls = auspice["controls"],
           geo_attributes = params["geo_inference"],
           color_options=auspice["color_options"],
           panels=auspice["panels"])
