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


region_cmap = [
    ["africa",          "#5097BA"],
    ["south_america",   "#60AA9E"],
    ["west_asia",       "#75B681"],
    ["oceania",         "#8EBC66"],
    ["europe",          "#AABD52"],
    ["japan_korea",     "#C4B945"],
    ["north_america",   "#D9AD3D"],
    ["southeast_asia",  "#E59637"],
    ["south_asia",      "#E67030"],
    ["china",           "#DF4327"]]

color_options = {
    "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete", "color_map": region_cmap},
    "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
}

date_range = {'date_min': '1920-01-01', 'date_max': '2017-06-01'}
geo_attributes = ['region']
attribute_nesting = {'geographic location':geo_attributes, 'authors':['authors']}
panels = ['tree', 'map', 'entropy']


def select_serotype(infile, path, serotype):
    '''
    Takes any-serotype fasta file, path to save output, and desired serotype.
    Writes appropriate subset of sequences to dengue_serotype.fasta
    Returns path to output file as string.
    '''
    sequences = [ i for i in SeqIO.parse(infile, 'fasta') if i.description.split('/')[0] == 'DENV%s'%serotype ]
    SeqIO.write(sequences, path+'dengue_%s.fasta'%serotype, 'fasta')
    return path+'dengue_%s.fasta'%serotype

class dengue_process(process):
    def __init__(self, **kwargs):
        super(process, self).__init__()

        self.serotype = kwargs['serotype']
        if self.serotype == 'any': # For any-serotype build, use dengue 4 outgroup and look for files like dengue.fasta
            self.lineage = 'dengue'
            self.reference_fname = './dengue/metadata/dengue_%s_outgroup.gb'%'4'
            newest_sequence_file = sorted(glob('../fauna/data/%s.fasta'%self.lineage), key=lambda f: os.path.getmtime(f))[-1]
        else:
            self.lineage = 'dengue_%s'%self.serotype # For serotype-specific build, use the corresponding outgroup
            self.reference_fname = './dengue/metadata/dengue_%s_outgroup.gb'%self.serotype
            try: # Look for a serotype-specific fasta
                newest_sequence_file = sorted(glob('../fauna/data/%s*.fasta'%self.lineage), key=lambda f: os.path.getmtime(f))[-1]
            except: # If it doesn't exist, try to pull serotype-specific sequences out of the any-serotype fasta (warn the user of this behavior)
                newest_sequence_file = select_serotype('../fauna/data/dengue.fasta', '../fauna/data/', self.serotype)
                print('WARNING: Did not find serotype-specific fasta file.\nPulled sequences with serotype %s from any-serotype fasta file %s\nWrote these to file %s'%(self.serotype, '../fauna/data/dengue.fasta', newest_sequence_file))

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
                       method='SLSQP',
                       lat_long_fname='../fauna/source-data/geo_lat_long.tsv')

        if params.load:
            self.dengue.load()
        else:
            self.fasta_fields = {0:'strain', 1:'accession', 2:'date', 3:'region', 4:'country',
                            5:'division', 6: 'location', 7: 'authors', 8: 'url'}
            self.dengue.load_sequences(fields=self.fasta_fields)
            self.dengue.seqs.filter(lambda s: len(s.seq)>=5000)
            self.dropped_strains = ['DENV1/VIETNAM/BIDV992/2006', 'DENV1/FRANCE/00475/2008', 'DENV1/VIETNAM/BIDV3990/2008', 'DENV2/HAITI/DENGUEVIRUS2HOMOSAPIENS1/2016'] # probable recombinants
            self.dengue.seqs.filter(lambda s: s.id not in self.dropped_strains)
            self.dengue.seqs.filter(lambda s: s.attributes['region'] not in ['', '?'])
            self.dengue.seqs.subsample(category = lambda x:(x.attributes['region'],
                                                     x.attributes['date'].year,
                                                     x.attributes['date'].month), threshold=params.viruses_per_month)
            for s in self.dengue.seqs.seqs.values():
                s.attributes['url'] = 'https://www.ncbi.nlm.nih.gov/nuccore/%s'%s.attributes['accession']
            self.dengue.align()
            self.dengue.build_tree()
            self.dengue.dump()
        self.dengue.clock_filter(n_iqd=3, plot=True)
        self.dengue.annotate_tree(Tc=False, timetree=True, reroot='best')

        for level in geo_attributes:
            self.dengue.tree.geo_inference(level)

        self.date_range = {'date_min': self.dengue.tree.getDateMin(), 'date_max': self.dengue.tree.getDateMax()}
        self.dengue.export(controls = attribute_nesting, geo_attributes = geo_attributes, date_range = self.date_range, color_options=color_options, panels=panels)


if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-v', '--viruses_per_month', type = int, default = 3, help='number of viruses sampled per month')
    parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
    parser.add_argument('--load', action='store_true', help = 'recover from file')
    parser.add_argument('-s', '--serotype', type = str, choices=['1', '2', '3', '4', 'any', 'all'], default='any', help = 'which serotype of dengue to build trees for; any = include all serotypes in one build; all = run all 5 builds')
    params = parser.parse_args()

    if params.serotype != 'all': # Run single build and quit.
        output = dengue_process(**params.__dict__)
    else: # Run all 5 builds. Call explicitly rather than loop to allow for ipython debugging of each build.
        setattr(params, 'serotype', '1')
        output1 = dengue_process(**params.__dict__)
        setattr(params, 'serotype', '2')
        output2 = dengue_process(**params.__dict__)
        setattr(params, 'serotype', '3')
        output3 = dengue_process(**params.__dict__)
        setattr(params, 'serotype', '4')
        output4 = dengue_process(**params.__dict__)
        setattr(params, 'serotype', 'any')
        outputany = dengue_process(**params.__dict__)
