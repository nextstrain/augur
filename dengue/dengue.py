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


class dengue_process(process):
    def __init__(self, **kwargs):
        super(process, self).__init__()

        self.serotype = kwargs['serotype']
        if self.serotype == 'all':
            self.lineage = 'dengue_virus'
            self.reference_fname = './dengue/metadata/dengue_virus_%s_outgroup.gb'%'3'
        else:
            self.lineage = 'dengue_virus_%s'%self.serotype
            self.reference_fname = './dengue/metadata/dengue_virus_%s_outgroup.gb'%self.serotype

        newest_sequence_file = sorted(glob('../fauna/data/%s*.fasta'%self.lineage))[-1]
        self.input_data_path = newest_sequence_file.split('.fasta')[0]
        self.sequence_fname = newest_sequence_file
        self.store_data_path = 'store/'+self.lineage + '_'
        self.build_data_path = 'build/'+self.lineage + '_'
        self.proteins = ['C', 'M', 'E', 'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', '2K', 'NS4B', 'NS5']
        self.load_reference(self.reference_fname, self.proteins)

        self.dengue = process(input_data_path = self.input_data_path,
                       store_data_path = self.store_data_path,
                       build_data_path = self.build_data_path,
                       proteins=self.proteins,
                       method='SLSQP',
                       reference_seq = self.reference_seq)

        if params.load:
            self.dengue.load()
        else:
            self.fasta_fields = {0:'strain', 1:'accession', 2:'date', 3:'region', 4:'country',
                            5:'division', 6: 'location'}
            self.dengue.load_sequences(fields=self.fasta_fields)


            self.dengue.seqs.filter(lambda s: s.attributes['date']>=datetime(2012,1,1).date() and
                                       s.attributes['date']< datetime(2017,1,1).date())
            self.dengue.seqs.filter(lambda s: len(s.seq)>=2000)
            self.dropped_strains = []
            self.dengue.seqs.filter(lambda s: s.id not in self.dropped_strains)
            self.dengue.seqs.subsample(category = lambda x:(x.attributes['region'],
                                                     x.attributes['date'].year,
                                                     x.attributes['date'].month), threshold=1000)

            self.dengue.align()
            self.dengue.build_tree()

        self.dengue.clock_filter(n_iqd=3, plot=True)
        self.dengue.annotate_tree(Tc=0.005, timetree=True, reroot='best')
        self.dengue.tree.geo_inference('region')
        self.dengue.tree.geo_inference('country')
        self.dengue.tree.geo_inference('division')
        self.dengue.export(controls = attribute_nesting)


    def load_reference(self, reference_file, protein_list):
        from Bio import SeqIO
        from Bio.SeqFeature import FeatureLocation
        refseq_names = {'1': 'DENV1/NAURUISLAND/REFERENCE/1997', '2': 'DENV2/THAILAND/REFERENCE/1964', '3': 'DENV3/SRI_LANKA/REFERENCE/2000', '4': 'DENV4/NA/REFERENCE/2003'}
        self.reference_seq = SeqIO.read(reference_file, 'genbank')
        self.reference_seq.id = refseq_names[self.serotype]
        self.reference_seq.name = refseq_names[self.serotype]
        self.genome_annotation = self.reference_seq.features
        if protein_list != []:
            # grab annotation from genbank
            self.proteins = {}
            for f in self.genome_annotation:
                if 'product' in f.qualifiers:
                    mat_peptide = f.qualifiers['product'][0].split()[-1]
                    location = FeatureLocation(start=f.location.start, end=f.location.end, strand=1)
                    if mat_peptide in protein_list and mat_peptide not in self.proteins: # keep the larger entry (e.g., keep prM instead of just pr and just M)
                        self.proteins[mat_peptide] = location
        else:
            self.proteins = {}

if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-v', '--viruses_per_month', type = int, default = 100, help='number of viruses sampled per month')
    parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
    parser.add_argument('--load', action='store_true', help = 'recover from file')
    parser.add_argument('-serotype', type = str, choices=['1', '2', '3', '4', 'all'], help = 'which serotype of dengue to build trees for')
    params = parser.parse_args()
    attribute_nesting = {'geographic location':['region', 'country', 'location']}
    dengue_process(**params.__dict__)
