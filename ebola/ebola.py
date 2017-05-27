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

attribute_nesting = {'geographic location':['country', 'division'], 'authors':['authors']}
date_range = {'date_min': '2013-01-01', 'date_max': '2017-01-01'}
geo_attributes = ['country', 'division']


country_cmap = [
    ["guinea",            "#86BB6D"],
    ["sierra_leone",      "#4041C7"],
    ["liberia",           "#E4682F"]
]
division_cmap = [
    # generic country colors
    ["guinea",            "#86BB6D"],
    ["sierra_leone",      "#4041C7"],
    ["liberia",           "#E4682F"],
    # Guinean prefectures
    ["boffa",             "#58A2AB"],
    ["boke",              "#56A0AE"],
    ["fria",              "#5BA7A6"],
    ["gaoual",            "#62AB9C"],
    ["koundara",          "#64AC99"],
    ["conakry",           "#60AA9F"],
    ["dabola",            "#84BA70"],
    ["dinguiraye",        "#81B973"],
    ["faranah",           "#86BB6D"],
    ["kissidougou",       "#95BD61"],
    ["kankan",            "#92BC63"],
    ["kerouane",          "#9EBE5A"],
    ["kouroussa",         "#89BB6B"],
    ["mandiana",          "#8FBC66"],
    ["siguiri",           "#8CBB68"],
    ["coyah",             "#66AE95"],
    ["dubreka",           "#5DA8A2"],
    ["forecariah",        "#68AF92"],
    ["kindia",            "#6AB18F"],
    ["telimele",          "#5AA4A8"],
    ["koubia",            "#7EB976"],
    ["labe",              "#74B582"],
    ["lelouma",           "#6FB388"],
    ["mali",              "#71B485"],
    ["tougue",            "#7CB879"],
    ["dalaba",            "#77B67F"],
    ["mamou",             "#79B77C"],
    ["pita",              "#6CB28C"],
    ["beyla",             "#A1BE58"],
    ["gueckedou",         "#98BD5E"],
    ["lola",              "#AABD52"],
    ["macenta",           "#9BBE5C"],
    ["nzerekore",         "#A4BE56"],
    ["yamou",             "#A7BE54"],
    # Sierra Leone districts
    ["kailahun",          "#4C90C0"],
    ["kenema",            "#447ECC"],
    ["kono",              "#4272CE"],
    ["bombali",           "#691D93"],
    ["kambia",            "#781C86"],
    ["koinadugu",         "#4041C7"],
    ["port_loko",         "#691D93"],
    ["tonkolili",         "#482BB6"],
    ["bo",                "#4066CF"],
    ["bonthe",            "#4041C7"],
    ["moyamba",           "#4041C7"],
    ["pujehun",           "#4887C6"],
    ["western_rural",     "#482BB6"],
    ["western_urban",     "#691D93"],
    ["western_area",      "#5824A4"],
    # Liberian counties
    ["nimba",             "#E67D33"],
    ["rivercess",         "#E68735"],
    ["river_gee",         "#E29E39"],
    ["sinoe",             "#E49838"],
    ["bomi",              "#DF4528"],
    ["bong",              "#E4682F"],
    ["gbarpolu",           "#DC2D24"],
    ["grand_cape_mount",  "#DB2122"],
    ["grand_bassa",       "#E67431"],
    ["grand_gedeh",       "#E69036"],
    ["grand_kru",         "#DFA53B"],
    ["lofa",              "#DE3926"],
    ["margibi",           "#E35C2C"],
    ["maryland",          "#DCAB3C"],
    ["montserrado",       "#E1512A"]
]

color_options = {
    "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete", "color_map":country_cmap},
    "division":{"key":"division", "legendTitle":"Division", "menuItem":"division", "type":"discrete", "color_map":division_cmap},
    "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
    "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"}
}
panels = ['tree', 'map', 'entropy']

if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-v', '--viruses_per_month', type = int, default = 100, help='number of viruses sampled per month')
    parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
    parser.add_argument('--load', action='store_true', help = 'recover from file')
    params = parser.parse_args()

    lineage = 'ebola'
    input_data_path = '../fauna/data/'+lineage
    store_data_path = 'store/'+lineage + '_'
    build_data_path = 'build/'+lineage + '_'

    ebola = process(input_data_path = input_data_path, store_data_path = store_data_path, build_data_path = build_data_path,
                   reference='ebola/metadata/ebola_outgroup.gb',
                   proteins=['NP', 'VP35', 'VP40', 'GP', 'sGP', 'VP30', 'VP24', 'L'],
                   method='SLSQP')

    if params.load==False:
        fasta_fields = {0:'strain', 2:'accession', 3:'date', 4:'region', 5:'country', 6:'division', 8:'db', 10:'authors', 11:'url'}
        ebola.load_sequences(fields=fasta_fields)
        ebola.seqs.filter(lambda s: s.attributes['date']>=datetime(2012,1,1).date() and
                                   s.attributes['date']< datetime(2017,1,1).date())
        dropped_strains = [
            "EM_074335", "EM_074462" # lack country metadata
        ]
        ebola.seqs.filter(lambda s: s.id not in dropped_strains)
        forced_strains = [
            "EM_076610" # flare-up index case
        ]
        ebola.seqs.subsample(category = lambda x:(x.attributes['region'], x.attributes['date'].year, x.attributes['date'].month),
            threshold=params.viruses_per_month, priority = lambda x:x.id in forced_strains)

        ebola.align()
        ebola.build_tree()
        ebola.dump()
    else:
        ebola.load()

    ebola.clock_filter(n_iqd=10, plot=True)
    ebola.annotate_tree(Tc="skyline", timetree=True, reroot='best', resolve_polytomies=True,
                        n_points=20, stiffness=3.0)
    for geo_attr in geo_attributes:
        ebola.tree.geo_inference(geo_attr)
    date_range['date_min'] = ebola.tree.getDateMin()
    date_range['date_max'] = ebola.tree.getDateMax()
    ebola.export(controls = attribute_nesting, date_range = date_range, geo_attributes = geo_attributes,
                color_options=color_options, panels=panels, defaults = {'mapTriplicate': False})

    # plot an approximate skyline
    plot_skyline = False
    if plot_skyline: # plot an approximate skyline
        from matplotlib import pyplot as plt
        T = ebola.tree.tt
        plt.figure()
        skyline, confidence = T.merger_model.skyline_inferred(gen = 50, confidence=2.0)
        plt.fill_between(skyline.x, confidence[0], confidence[1], color=(0.8, 0.8, 0.8))
        plt.plot(skyline.x, skyline.y)
        plt.yscale('log')
        plt.ylabel('effective population size')
        plt.xlabel('year')
        plt.ticklabel_format(axis='x',useOffset=False)
        plt.savefig('ebola_skyline.png')
