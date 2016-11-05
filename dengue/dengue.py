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

attribute_nesting = {'geographic location':['region', 'country', 'location'],}

class dengue_process(process):
    def __init__(self, **kwargs):
        self.sequence_fname = self.input_data_path+'/dengue_virus_%d.fasta'%self.serotype

    def load_reference(self, reference_file):
        from Bio import SeqIO
        from Bio.SeqFeature import FeatureLocation
        self.reference_seq = SeqIO.read(reference_file, 'genbank')
        self.reference_seq.id = self.reference_seq.name
        self.genome_annotation = self.reference_seq.features
        if "proteins" in self.kwargs:
            # grap annotation from genbank
            protein_list = self.kwargs['proteins']
            self.proteins = {f.qualifiers['gene'][0]:FeatureLocation(start=f.location.start, end=f.location.end, strand=1)
                            for f in self.genome_annotation
                                if 'gene' in f.qualifiers
                                    and f.qualifiers['gene'][0] in protein_list}
        else:
            self.proteins = {}



if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-v', '--viruses_per_month', type = int, default = 100, help='number of viruses sampled per month')
    parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
    parser.add_argument('--load', action='store_true', help = 'recover from file')
    parser.add_argument('-serotype', type = str, choices=['1', '2', '3', '4'], help = 'which serotype of dengue to build trees for')
    params = parser.parse_args()

    lineage = 'dengue'
    input_data_path = '../fauna/data/'+lineage ## Should add an assertion here that these directories exist.
    store_data_path = 'store/'+lineage + '_'
    build_data_path = 'build/'+lineage + '_'

    dengue = process(input_data_path = input_data_path,
                   store_data_path = store_data_path,
                   build_data_path = build_data_path,
                   reference='./dengue/metadata/dengue_virus_%s_outgroup.gb'%params.serotype,
                   proteins=['C', 'prM', 'E', 'NS1', '2A', '2B',
                             'NS3', '4A', '4B', 'NS5'],
                   method='SLSQP')

    if params.load:
        dengue.load()
    else:
        fasta_fields = {0:'strain', 2:'accession', 3:'date', 4:'region', 5:'country',
                        6:'division', 8:'db', 10:'authors', 11:'latitude', 12:'longitude'}
        dengue.load_sequences(fields=fasta_fields)
        dengue.seqs.filter(lambda s: s.attributes['date']>=datetime(2012,1,1).date() and
                                   s.attributes['date']< datetime(2017,1,1).date())
        dengue.seqs.filter(lambda s: len(s.seq)>=2000)
        dropped_strains = []
        dengue.seqs.filter(lambda s: s.id not in dropped_strains)
        dengue.seqs.subsample(category = lambda x:(x.attributes['region'],
                                                 x.attributes['date'].year,
                                                 x.attributes['date'].month), threshold=1000)

        dengue.align()
        dengue.build_tree()

    dengue.clock_filter(n_iqd=3, plot=True)
    dengue.annotate_tree(Tc=0.005, timetree=True, reroot='best')
    dengue.tree.geo_inference('region')
    dengue.tree.geo_inference('country')
    dengue.tree.geo_inference('division')
    dengue.export(controls = attribute_nesting)
