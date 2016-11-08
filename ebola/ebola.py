from __future__ import division, print_function
from collections import defaultdict
import sys
sys.path.append('')  # need to import from base
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
geo_attributes = ['country', 'division']

color_options = {
    "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
    "division":{"key":"division", "legendTitle":"District", "menuItem":"District", "type":"discrete"},
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
        fasta_fields = {0:'strain', 2:'accession', 3:'date', 4:'region', 5:'country', 6:'division', 8:'db', 10:'authors'}
        ebola.load_sequences(fields=fasta_fields)
        ebola.seqs.filter(lambda s: s.attributes['date']>=datetime(2012,1,1).date() and
                                   s.attributes['date']< datetime(2017,1,1).date())
        forced_strains = [
            "EM_076610" # flare-up index case
        ]
        ebola.seqs.subsample(category = lambda x:(x.attributes['region'], x.attributes['date'].year, x.attributes['date'].month),
            threshold=params.viruses_per_month, forced_strains = forced_strains)

        ebola.align()
        ebola.build_tree()
        ebola.dump()
    else:
        ebola.load()

    ebola.clock_filter(n_iqd=10, plot=True)
    ebola.annotate_tree(Tc=0.0005, timetree=True, reroot='best', resolve_polytomies=True)
    for geo_attr in geo_attributes:
        ebola.tree.geo_inference(geo_attr)
    ebola.export(controls = attribute_nesting, geo_attributes = geo_attributes,
                color_options=color_options, panels=panels)
