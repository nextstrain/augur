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

attribute_nesting = {'geographic location':['region', 'country', 'city']}

if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-v', '--viruses_per_month', type = int, default = 20, help='number of viruses sampled per month')
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
    #    dropped_strains = []
    #    ebola.seqs.filter(lambda s: s.id not in dropped_strains)
        ebola.seqs.subsample(category = lambda x:(x.attributes['region'],
                                                 x.attributes['date'].year,
                                                 x.attributes['date'].month), threshold=par.viruses_per_month)

        ebola.align()
        ebola.build_tree()
        ebola.dump()
    else:
        ebola.load()

    ebola.clock_filter(n_iqd=10, plot=True)
    ebola.annotate_tree(Tc=0.0005, timetree=True, reroot='best', resolve_polytomies=True)
    #ebola.tree.geo_inference('region')
    ebola.tree.geo_inference('country')
    ebola.tree.geo_inference('division')
    ebola.export(controls = attribute_nesting)

