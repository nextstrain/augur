from __future__ import division, print_function
from collections import defaultdict
from base.io_util import make_dir, remove_dir, tree_to_json, write_json, myopen
from base.sequences import sequence_set, num_date
from base.tree import tree
from base.process import process
from base.frequencies import alignment_frequencies, tree_frequencies
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import numpy as np
from datetime import datetime

if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-v', '--viruses_per_month', type = int, default = 100, help='number of viruses sampled per month')
    parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
    parser.add_argument('--load', action='store_true', help = 'recover from file')

    lineage = 'zika'
    input_data_path = '../nextstrain-db/data/'+lineage
    store_data_path = 'store/'+lineage + '_'
    build_data_path = 'build/'+lineage + '_'

    zika = process(input_data_path = input_data_path, store_data_path = store_data_path, build_data_path = build_data_path,
                   reference='zika/metadata/zika_outgroup.gb',
                   proteins=['CA', 'PRO', 'MP', 'ENV', 'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', 'NS4B'],
                   method='SLSQP')

    fasta_fields = {0:'strain', 2:'isolate_id', 3:'date', 4:'region',
                    5:'country', 7:"city"}
    zika.load_sequences(fields=fasta_fields)
    zika.seqs.filter(lambda s: s.attributes['date']>=datetime(2012,1,1).date() and
                               s.attributes['date']< datetime(2017,1,1).date())
    zika.seqs.subsample(category = lambda x:(x.attributes['region'],
                                             x.attributes['date'].year,
                                             x.attributes['date'].month), threshold=100)

    zika.align()
    zika.build_tree()
    zika.annotate_tree(Tc=0.005, timetree=True, reroot='clock_filter', n_iqd=3, plot=True)
    zika.tree.geo_inference('region')
    zika.export()
