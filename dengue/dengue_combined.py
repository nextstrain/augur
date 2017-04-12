from __future__ import division, print_function
import sys
sys.path.append('')  # need to import from base
from base.io_util import make_dir, remove_dir, tree_to_json, write_json, myopen
from base.sequences import sequence_set, num_date
from base.tree import tree
from base.process import process
from base.frequencies import alignment_frequencies, tree_frequencies
from Bio import SeqIO
# from Bio.SeqFeature import FeatureLocation
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

regions = ['africa', 'south_asia', 'europe', 'china', 'north_america',
           'china', 'south_america', 'japan_korea', 'oceania', 'southeast_asia']

region_groups = {'NA':'north_america',
                 'AS':['china', 'japan_korea', 'south_asia', 'southeast_asia'],
                 'OC':'oceania', 'EU':'europe'}

color_options = {
    "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete", "color_map": region_cmap},
    "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
    "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"},
    # "titers": {"key": "titers", "legendTitle": "Titers", "menuItem": "titers", "type":"continuous"}
}

# date_range = {'date_min': '1920-01-01', 'date_max': '2017-06-01'}
attribute_nesting = {'geographic location':'region', 'authors':['authors']}
panels = ['tree', 'map', 'entropy']#, 'frequencies']

pivots = np.arange(1920.,2017.)

genotypes = {
                        "dengue_1":{
                        'I': [('', 961, 'G'),('', 965, 'G'),('', 995, 'T'),('', 1001, 'A'),('', 1016, 'A')],
                        'II': [('', 995, 'A'),('', 1100, 'A'),('', 1103, 'T'),('', 1109, 'T'),('', 1113, 'T')],
                        "V": [('',974,'T'), ('', 1049,'T'),('',1175,'G'),('',1091,'G')],
                        "IV": [('', 944,'D'), ('', 992,'T'),('', 1034,'T'),('', 1047,'G'),('', 1181,'G')]},

                        "dengue_2": {
                        'ASIAN/AMERICAN': [('', 1181, 'T'),('', 974, 'A'),('', 1025, 'T'),('', 1103, 'T')],
                        'AMERICAN': [('', 965, 'G'),('', 983, 'G'),('', 992, 'G'),('', 995, 'T')],
                        'ASIAN-I': [('', 974, 'A'),('', 995, 'T'),('', 1067, 'T'),('', 1082, 'G')],
                        'COSMOPOLITAN': [('', 947, 'T'),('', 1150, 'C'),('', 1253, 'C'),('', 1277, 'C')],
                        'SYLVATIC': [('', 977, 'G'),('', 1100, 'A'),('', 1103, 'A')]},

                        "dengue_3": {
                        'I': [('', 980, 'C'),('', 1034, 'T'),('', 1073, 'C'),('', 1118, 'T')],
                        'II': [('', 983, 'G'),('', 1004, 'T'),('', 1064, 'C'),('', 1113, 'T'),('', 1232, 'T')],
                        'III': [('', 1022, 'G'),('', 1127, 'G'),('', 1157, 'A')],
                        'IV': [('', 956, 'G'),('', 980, 'C'),('', 986, 'G'),('', 1001, 'C')]},

                        "dengue_4": {
                        'I': [('', 974, 'A'),('', 1004, 'T'),('', 1130, 'T'),('', 1181, 'C'),('', 1199, 'T')],
                        'II': [('', 1075, 'C'),('', 1223, 'T'),('', 1322, 'T'),('', 1433, 'T')],
                        'SYLVATIC': [('', 956, 'G'),('', 962, 'T'),('', 1019, 'C')]},

                        "all": {}
                }
# for i in ['dengue_1', 'dengue_2', 'dengue_3', 'dengue_4']:
#     genotypes['all'].update(genotypes[i])

def select_serotype(infile, path, serotype):
    '''
    Takes all-serotype fasta file, path to save output, and desired serotype.
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
        if self.serotype == 'all': # For all-serotype build, use dengue 4 outgroup and look for files like dengue.fasta
            self.lineage = 'dengue_all'
            self.reference_fname = './dengue/metadata/dengue_%s_outgroup.gb'%'4'
            newest_sequence_file = sorted(glob('../fauna/data/dengue.fasta'), key=lambda f: os.path.getmtime(f))[-1]
        else:
            self.lineage = 'dengue_denv%s'%self.serotype # For serotype-specific build, use the corresponding outgroup
            self.reference_fname = './dengue/metadata/dengue_%s_outgroup.gb'%self.serotype
            try: # Look for a serotype-specific fasta
                newest_sequence_file = sorted(glob('../fauna/data/%s*.fasta'%self.lineage), key=lambda f: os.path.getmtime(f))[-1]
            except: # If it doesn't exist, try to pull serotype-specific sequences out of the all-serotype fasta (warn the user of this behavior)
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
                       method='SLSQP',
                       lat_long_fname='../fauna/source-data/geo_lat_long.tsv')

        if params.load:
            self.dengue.load()
        else:
            self.fasta_fields = {0:'strain', 1:'accession', 2:'date', 3:'region', #4:'country',
                            5:'division', 6: 'location', 7: 'authors', 8: 'url'}
            self.dengue.load_sequences(fields=self.fasta_fields)
            self.dengue.seqs.filter(lambda s: len(s.seq)>=5000)
            self.dropped_strains = ['DENV1/VIETNAM/BIDV992/2006', 'DENV1/FRANCE/00475/2008', 'DENV1/VIETNAM/BIDV3990/2008', 'DENV2/HAITI/DENGUEVIRUS2HOMOSAPIENS1/2016', # Probable recombinants
            'DENV2/AUSTRALIA/QML22/2015'] # Suspiciously far diverged
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

        nodes = [i for i in self.dengue.tree.tree.find_clades() ]
        print(nodes[110].attr)

        self.dengue.annotate_tree(Tc=False, timetree=True, reroot='best')
        nodes = [i for i in self.dengue.tree.tree.find_clades() ]
        print(nodes[110].attr)

        # self.dengue.estimate_mutation_frequencies(region="global", pivots=pivots, min_freq=.01)
        # for region in region_groups.iteritems():
        #     self.dengue.estimate_mutation_frequencies(region=region, min_freq=.05)

        # self.dengue.estimate_tree_frequencies()
        # for region in regions:
        #     self.dengue.estimate_tree_frequencies(region=region)

        # self.dengue.matchClades(genotypes[self.lineage])

        for level in geo_attributes:
            self.dengue.tree.geo_inference(level)

        nodes = [i for i in self.dengue.tree.tree.find_clades() ]
        print(nodes[110].attr)

        self.date_range = {'date_min': self.dengue.tree.getDateMin(), 'date_max': self.dengue.tree.getDateMax()}
        self.dengue.export(controls = attribute_nesting, geo_attributes = geo_attributes, date_range = self.date_range, color_options=color_options, panels=panels)


if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-v', '--viruses_per_month', type = int, default = 3, help='number of viruses sampled per month')
    parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
    parser.add_argument('-s', '--serotype', type = str, choices=['1', '2', '3', '4', 'all'], default='all', help = 'which serotype of dengue to build trees for; all = include all serotypes in one build')
    parser.add_argument('--load', action='store_true', help = 'recover from file')
    parser.add_argument('--process_multiple', action='store_true', help = 'run all serotype builds')
    params = parser.parse_args()

    if params.process_multiple: # Run all 5 builds. Call explicitly rather than loop to allow for ipython debugging of each build.
        setattr(params, 'serotype', '1')
        output1 = dengue_process(**params.__dict__)
        setattr(params, 'serotype', '2')
        output2 = dengue_process(**params.__dict__)
        setattr(params, 'serotype', '3')
        output3 = dengue_process(**params.__dict__)
        setattr(params, 'serotype', '4')
        output4 = dengue_process(**params.__dict__)
        setattr(params, 'serotype', 'all')
        outputall = dengue_process(**params.__dict__)
    else: # Run single build and quit.
        output = dengue_process(**params.__dict__)
