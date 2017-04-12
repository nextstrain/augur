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


color_options = {
    "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
}

date_range = {'date_min': '1920-01-01', 'date_max': '2017-06-01'}
panels = ['tree']

class dengue_process(process):
    def __init__(self, **kwargs):
        super(process, self).__init__()

        self.lineage = 'dengue'
        self.reference_fname = './dengue/metadata/dengue_%s_outgroup.gb'%'4'
        self.sequence_fname = '../fauna/data/genotype_codonAln_E.fasta'
        self.input_data_path = '../fauna/data/genotype_codonAln_E'
        self.store_data_path = 'store/genotypes_'
        self.build_data_path = 'build/genotypes_'
        self.proteins = ['E']

        self.dengue = process(input_data_path = self.input_data_path,
                       store_data_path = self.store_data_path,
                       build_data_path = self.build_data_path,
                       proteins=self.proteins,
                       reference=self.reference_fname,
                       method='SLSQP')

        self.fasta_fields = {0:'strain', 1:'accession'}
        self.dengue.load_sequences(fields=self.fasta_fields, prune=False)
        for seq in self.dengue.seqs.all_seqs.values():
            seq.attributes['date'] = seq.id.split('/')[3].split('|')[-1]
        self.dengue.seqs.parse_date(["%Y"], prune=False)
        self.dengue.seqs.seqs = self.dengue.seqs.all_seqs #subsample(threshold=45)
        self.dengue.align()
        self.dengue.build_tree()
        self.dengue.dump()
        self.dengue.clock_filter(n_iqd=3, plot=True)
        self.dengue.annotate_tree(Tc=False, timetree=True, reroot='best')

        self.date_range = {'date_min': self.dengue.tree.getDateMin(), 'date_max': self.dengue.tree.getDateMax()}
        self.dengue.export( date_range = self.date_range, color_options=color_options, panels=panels)


if __name__=="__main__":
    output = dengue_process()
