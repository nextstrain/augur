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
import os
from glob import glob

def select_serotype(infile, path, serotype):
    firstfield = 'DENV%s'%serotype
    sequences = [ i for i in SeqIO.parse(infile, 'fasta') if i.description.split('/')[0] == firstfield ]
    SeqIO.write(sequences, path+'dengue_%s.fasta'%serotype, 'fasta')
    return path+'dengue_%s.fasta'%serotype

class dengue_process(process):
    def __init__(self, **kwargs):
        super(process, self).__init__()

        self.serotype = kwargs['serotype']
        if self.serotype == 'all': # For any-serotype build, use dengue 3 outgroup and look for files like dengue.fasta
            self.lineage = 'dengue'
            self.reference_fname = './dengue/metadata/dengue_%s_outgroup.gb'%'3'
            newest_sequence_file = sorted(glob('../fauna/data/%s.fasta'%self.lineage), key=lambda f: os.path.getmtime(f))[-1]
        else:
            self.lineage = 'dengue_%s'%self.serotype # For serotype-specific build, use the corresponding outgroup
            self.reference_fname = './dengue/metadata/dengue_%s_outgroup.gb'%self.serotype
            try: # Look for a serotype-specific fasta
                newest_sequence_file = sorted(glob('../fauna/data/%s*.fasta'%self.lineage), key=lambda f: os.path.getmtime(f))[-1]
            except: # If it doesn't exist, try to pull serotype-specific sequences out of the any-serotype fasta (warn the user of this behavior)
                newest_sequence_file = select_serotype('../fauna/data/dengue.fasta', '../fauna/data/', self.serotype)
                print('WARNING: Did not find serotype-specific fasta file.\nPulled sequences with serotype %s from all-serotype fasta file %s\nWrote these to file %s'%(self.serotype, '../fauna/data/dengue.fasta', newest_sequence_file))

        self.input_data_path = newest_sequence_file.split('.fasta')[0]
        self.sequence_fname = newest_sequence_file
        self.store_data_path = 'store/'+self.lineage + '_'
        self.build_data_path = 'build/'+self.lineage + '_'
        self.proteins = ['C', 'M', 'E', 'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', '2K', 'NS4B', 'NS5']

        self.dengue = process(input_data_path = self.input_data_path,
                       store_data_path = self.store_data_path,
                       build_data_path = self.build_data_path,
                       proteins=self.proteins,
                       reference=self.reference_fname,
                       method='SLSQP')

        if params.load:
            self.dengue.load()
        else:
            self.fasta_fields = {0:'strain', 1:'accession', 2:'date', 3:'region', 4:'country',
                            5:'division', 6: 'location'}
            self.dengue.load_sequences(fields=self.fasta_fields)
            self.dengue.seqs.filter(lambda s: len(s.seq)>=2000)
            self.dropped_strains = []
            self.dengue.seqs.filter(lambda s: s.id not in self.dropped_strains)
            self.dengue.seqs.subsample(category = lambda x:(x.attributes['region'],
                                                     x.attributes['date'].year,
                                                     x.attributes['date'].month), threshold=params.viruses_per_month)
            self.dengue.align()
            self.dengue.build_tree()

        self.dengue.clock_filter(n_iqd=3, plot=True)
        self.dengue.annotate_tree(Tc=0.005, timetree=True, reroot='best')
        self.dengue.tree.geo_inference('region')
        self.dengue.export(controls = attribute_nesting)

if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-v', '--viruses_per_month', type = int, default = 100, help='number of viruses sampled per month')
    parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
    parser.add_argument('--load', action='store_true', help = 'recover from file')
    parser.add_argument('-s', '--serotype', type = str, choices=['1', '2', '3', '4', 'all'], default='all', help = 'which serotype of dengue to build trees for')
    params = parser.parse_args()
    attribute_nesting = {'geographic location':['region', 'country', 'location']}
    output = dengue_process(**params.__dict__)
