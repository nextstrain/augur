from __future__ import division, print_function
from collections import defaultdict
import sys
sys.path.insert(0,'.')  # need to import from base
from base.io_util import make_dir, remove_dir, tree_to_json, write_json, myopen
from base.sequences import sequence_set, num_date
from base.tree import tree
from base.process import process
from base.frequencies import alignment_frequencies, tree_frequencies
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import numpy as np
from datetime import datetime

region_cmap = [
    ["europe",             "#4580CA"],
    # ["china",              "#3F63CF"],
    ["japan_korea",        "#4272CE"],
    ["southeast_asia",     "#571EA2"],
    # ["oceania",            "#4B26B1"],
    ["south_pacific",      "#4433BE"],
    ["south_america",      "#56A0AE"],
    ["central_america",    "#BDBB48"],
    ["caribbean",          "#E67C32"],
    ["north_america",      "#DC2F24"]
]

country_cmap = [

    # "#571EA2", "#4B26B1", "#4433BE", "#4042C7", "#3F52CD", "#3F63CF", "#4272CE", "#4580CA", "#4A8CC2", "#5097BA"
    # "#56A0AE", "#5DA8A3", "#66AE96", "#6EB389", "#79B77D", "#83BA70", "#8EBC66", "#9ABE5C", "#A6BE55", "#B2BD4D"
    # "#BDBB48", "#C8B944", "#D1B340", "#D9AD3D", "#DFA43B", "#E49938", "#E68B35", "#E67C32", "#E56A2F", "#E2562B", "#DF4327", "#DC2F24"

    # Old World: "#571EA2", "#4B26B1", "#4433BE", "#4042C7", "#3F52CD", "#3F63CF", "#4272CE", "#4580CA", "#4A8CC2", "#5097BA"
    # ["thailand",           "#571EA2"],
    #["vietnam",            "#"],
    ["singapore",          "#4B26B1"],
    ["french_polynesia",   "#4433BE"],
    ["american_samoa",	   "#4042C7"],
    ["fiji",	           "#3F52CD"],
    ["tonga",	           "#3F63CF"],
    ["italy",			   "#4272CE"],
    # ["china",			   "#4580CA"],
    ["japan",			   "#4A8CC2"],
    # South America: "#56A0AE", "#5DA8A3", "#66AE96", "#6EB389", "#79B77D", "#83BA70", "#8EBC66", "#9ABE5C"
    ["brazil",  		   "#56A0AE"],
    ["peru",               "#5DA8A3"],
    ["ecuador",  		   "#66AE96"],
    ["colombia",  		   "#6EB389"],
    ["french_guiana",  	   "#79B77D"],
    ["suriname",  		   "#83BA70"],
    ["venezuela",  		   "#8EBC66"],
    # Central and North America: "#A6BE55", "#B2BD4D", "#BDBB48", "#C8B944", "#D1B340", "#D9AD3D", "#DFA43B",
    # "#E49938", "#E68B35", "#E67C32", "#E56A2F", "#E2562B", "#DF4327", "#DC2F24"
    ["panama",             "#A6BE55"],
    ["nicaragua",          "#B2BD4D"],
    ["honduras",           "#BDBB48"],
    ["el_salvador",        "#C8B944"],
    ["guatemala",  		   "#D1B340"],
    ["mexico",             "#D9AD3D"],
    ["martinique",  	   "#DFA43B"],
    ["guadeloupe",         "#E49938"],
    ["usvi",               "#E68B35"],
    ["puerto_rico",  	   "#E67C32"],
    ["jamaica",            "#E56A2F"],
    ["dominican_republic", "#E2562B"],
    ["haiti",  			   "#DF4327"],
    ["usa",                "#DC2F24"]
]


attribute_nesting = {'geographic location':['region', 'country'], 'authors':['authors']}
date_range = {'date_min': '2013-06-01', 'date_max': '2017-06-01'}
geo_attributes = ['region', 'country']
panels = ['tree', 'map', 'entropy']
color_options = {
    "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete", "color_map": country_cmap},
    "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete", "color_map": region_cmap},
    "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
    "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"}
}

if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-v', '--viruses_per_month', type = int, default = 100, help='number of viruses sampled per month')
    parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
    parser.add_argument('--verbose', type = int, default = 2, help='throttle treetime output')
    parser.add_argument('--confidence', default = False, action="store_true", help='evaluate confidence intervals of internal node timing')
    parser.add_argument('--load', action='store_true', help = 'recover from file')
    params = parser.parse_args()

    lineage = 'zika'
    input_data_path = '../fauna/data/'+lineage
    store_data_path = 'store/'+lineage + '_'
    build_data_path = 'build/'+lineage + '_'

    zika = process(input_data_path = input_data_path,
                   store_data_path = store_data_path,
                   build_data_path = build_data_path,
                   reference='zika/metadata/zika_outgroup.gb',
                   lat_long_fname='../fauna/source-data/geo_lat_long.tsv',
                   proteins=['CA', 'PRO', 'MP', 'ENV', 'NS1', 'NS2A',
                             'NS2B', 'NS3', 'NS4A', 'NS4B', 'NS5'],
                   method='SLSQP', verbose=params.verbose)
    if params.load:
        zika.load()
    else:
        fasta_fields = {0:'strain', 2:'accession', 3:'date', 4:'region', 5:'country',
                        6:'division', 8:'db', 10:'authors', 11:'url'}
        zika.load_sequences(fields=fasta_fields)
        zika.seqs.filter(lambda s: s.attributes['date']>=datetime(2012,1,1).date() and
                                   s.attributes['date']< datetime(2018,1,1).date())
        zika.seqs.filter(lambda s: len(s.seq)>=2000)
        dropped_strains = [
            "THA/PLCal_ZV/2013", "PLCal_ZV", # true strains, too basal for analysis
            "ZF36_36S", # possible contamination
            "Dominican_Republic/2016/PD2", "GD01", "GDZ16001", "VEN/UF_2/2016", # true strains, but duplicates of other strains in dataset
            "Bahia04", "JAM/2016/WI_JM6", # excessive terminal branch length
            "THA/2014/SV0127_14", "ZK_YN001", "NIID123/2016", # true strains, too basal for analysis
            "ZKA_16_291", "ZKA_16_097" # singapore, too basal for analysis
        ]
        zika.seqs.filter(lambda s: s.id not in dropped_strains)
        zika.seqs.subsample(category = lambda x:(x.attributes['date'].year, x.attributes['date'].month),
            threshold=params.viruses_per_month)

        zika.align()
        zika.build_tree()
        zika.dump()
    zika.clock_filter(n_iqd=3, plot=True)
    zika.annotate_tree(Tc=0.02, timetree=True, reroot='best', confidence=params.confidence)
    zika.tree.geo_inference('region')
    zika.tree.geo_inference('country')
    date_range['date_min'] = zika.tree.getDateMin()
    date_range['date_max'] = zika.tree.getDateMax()
    zika.export(controls = attribute_nesting, date_range = date_range, geo_attributes = geo_attributes,
                color_options=color_options, panels=panels)

    plot_skyline = False
    if plot_skyline: # plot an approximate skyline
        from matplotlib import pyplot as plt
        T = zika.tree.tt
        plt.figure()
        skyline, confidence = T.merger_model.skyline_inferred(gen = 50, confidence=2.0)
        plt.fill_between(skyline.x, confidence[0], confidence[1], color=(0.8, 0.8, 0.8))
        plt.plot(skyline.x, skyline.y)
        plt.yscale('log')
        plt.ylabel('effective population size')
        plt.xlabel('year')
        plt.ticklabel_format(axis='x',useOffset=False)
        plt.savefig('zika_skyline.png')
