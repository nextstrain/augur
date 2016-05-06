from __future__ import division, print_function
import os, time, gzip
from collections import defaultdict
from nextstrain.io_util import make_dir, remove_dir, tree_to_json, write_json, myopen
from nextstrain.sequences import sequence_set, num_date
from nextstrain.tree import tree
from nextstrain.frequencies import alignment_frequencies, tree_frequencies
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import numpy as np
from datetime import datetime

def fix_name(name):
    tmp_name = name.replace('_', '').replace(' ', '').replace('\'','')\
                   .replace('(','').replace(')','').replace('H3N2','')\
                   .replace('Human','').replace('human','').replace('//','/')
    return tmp_name.strip().strip('_')

class zika_process(object):
    """process zika virus sequences in mutliple steps to allow visualization in browser
        * filtering and parsing of sequences
        * alignment
        * tree building
        * frequency estimation of clades and mutations
        * export as json
    """

    def __init__(self, fname = 'data/zika_full.fasta',
                 out_specs={'data_dir':'data/', 'prefix':'zika_', 'qualifier':''},
                 **kwargs):
        super(zika_process, self).__init__()
        self.fname = fname
        self.kwargs = kwargs
        self.out_specs = out_specs
        if 'outgroup' in kwargs:
            outgroup_file = kwargs['outgroup']
        else:
            outgroup_file = 'source_data/'+out_specs['prefix']+'outgroup2.gb'
        outliers = ['KU744693.1', 'KU866423.1', 'KU820897.1', 'KU940228','?']
        tmp_outgroup = SeqIO.read(outgroup_file, 'genbank')
        self.outgroup = tmp_outgroup.id.split('.')[0]
        genome_annotation = tmp_outgroup.features
        ref_seq = SeqIO.read(outgroup_file, 'genbank')
        self.proteins = {f.qualifiers['product'][0].split()[-1]:FeatureLocation(start=f.location.start, end=f.location.end, strand=1)
                for f in ref_seq.features if 'product' in f.qualifiers and f.type == 'mat_peptide'}

        self.time_interval = [datetime.strptime('2012-05-01', "%Y-%m-%d").date(),
                              datetime.strptime('2016-06-01', "%Y-%m-%d").date()]
        self.pivots = np.linspace(num_date(self.time_interval[0]),
                                  num_date(self.time_interval[1]),40)

        self.seqs = sequence_set(self.fname, reference=self.outgroup)
        print('found',len(self.seqs.raw_seqs), 'sequences')
        self.seqs.ungap()
        self.seqs.parse({0:'strain', 2:'isolate_id', 4:'region', 5:'country', 3:'date'}, strip='_')
        self.fix_strain_names()
        acc_to_strain = {s.attributes['isolate_id']:s.attributes['strain'] for s in self.seqs.raw_seqs.values()}
        strain_to_acc = {s.attributes['strain']:s.attributes['isolate_id'] for s in self.seqs.raw_seqs.values()}
        for seq in self.seqs.raw_seqs.values():
            seq.attributes['date'] = seq.attributes['date'].replace('XX', '01')
        self.seqs.raw_seqs[acc_to_strain[self.outgroup]].seq=tmp_outgroup.seq
        self.seqs.raw_seqs = {k:v for k,v in self.seqs.raw_seqs.iteritems() if k!=''}
        self.seqs.reference = self.seqs.raw_seqs[acc_to_strain[self.outgroup]]
        self.seqs.parse_date(["%Y-%m-%d", "%Y-%m", "%Y"], prune=True)
        print('after date parsing:',len(self.seqs.raw_seqs), 'sequences')
        self.seqs.raw_seqs = {k:s for k,s in self.seqs.raw_seqs.iteritems() if
                                        s.attributes['date']>=self.time_interval[0] and
                                        s.attributes['date']<self.time_interval[1] and
                                        strain_to_acc[k] not in outliers}
        self.filenames()


    def filenames(self):
        data_path = self.out_specs['data_dir']+self.out_specs['prefix']
        data_path += self.out_specs['qualifier']
        self.file_dumps = {}
        self.file_dumps['seqs'] = data_path+'sequences.pkl.gz'
        self.file_dumps['tree'] = data_path+'tree.newick'
        self.file_dumps['frequencies'] = data_path+'frequencies.pkl.gz'
        self.file_dumps['tree_frequencies'] = data_path+'tree_frequencies.pkl.gz'


    def fix_strain_names(self):
        new_raw_seqs = {}
        for desc, seq in self.seqs.raw_seqs.iteritems():
            new_name = fix_name(seq.attributes['strain'])
            seq.attributes['strain']=new_name
            seq.name=str(new_name)
            seq.id=str(new_name)
            new_raw_seqs[new_name]=seq
        self.seqs.raw_seqs = new_raw_seqs


    def dump(self):
        from cPickle import dump
        from Bio import Phylo
        for attr_name, fname in self.file_dumps.iteritems():
            if hasattr(self,attr_name):
                print("dumping",attr_name)
                if attr_name=='seqs': self.seqs.raw_seqs = None
                with myopen(fname, 'wb') as ofile:
                    if attr_name=='tree':
                        Phylo.write(self.tree.tree, ofile, 'newick')
                    else:
                        dump(getattr(self,attr_name), ofile, -1)

    def load(self):
        from cPickle import load
        for attr_name, fname in self.file_dumps.iteritems():
            if os.path.isfile(fname):
                with myopen(fname, 'r') as ifile:
                    if attr_name=='tree':
                        continue
                    else:
                        setattr(self, attr_name, load(ifile))
        fname = self.file_dumps['tree']
        if os.path.isfile(fname):
            self.build_tree(fname)

    def subsample(self):
        self.seqs.subsample(category = lambda x:(x.attributes['date'].year,
                                                 x.attributes['date'].month),
                            threshold=params.viruses_per_month)

    def align(self):
        self.seqs.align()
        self.seqs.strip_non_reference()
        #self.seqs.clock_filter(n_iqd=5, plot=True, max_gaps=0.05, root_seq=self.outgroup)
        for seq in self.seqs.aln:
            seq.seq = '-'*99 + seq.seq[99:]
        self.seqs.translate(proteins=self.proteins)

    def build_tree(self, infile=None):
        self.tree = tree(aln=self.seqs.aln, proteins = self.proteins)
        if infile is None:
            self.tree.build()
        else:
            self.tree.tt_from_file(infile)
        self.tree.timetree(Tc=0.005, infer_gtr=True)
        for node in self.tree.tt.tree.get_terminals():
            if hasattr(node, "bad_branch") and node.bad_branch:
                node.numdate = min(2016.15, node.numdate)
        self.tree.add_translations()
        self.tree.refine()
        self.tree.layout()


    def export(self, prefix='web/data/'):
        def process_freqs(freq):
            return [round(x,4) for x in freq]
        for n in self.tree.tree.find_clades():
            if hasattr(n,'muts'):
                n.nuc_muts= n.muts
        self.seqs.export_diversity(prefix+'entropy.json')
        self.tree.export(path=prefix, extra_attr = ["country", "region", "nuc_muts",
                                                    "aa_muts"])


if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-v', '--viruses_per_month', type = int, default = 100, help='number of viruses sampled per month')
    parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
    parser.add_argument('--load', action='store_true', help = 'recover from file')

    params = parser.parse_args()

    zika = zika_process(method='SLSQP', dtps=2.0, stiffness=20, inertia=0.9)
    if params.load:
        zika.load()
    else:
        zika.subsample()
        zika.align()
        zika.dump()
        zika.build_tree()
        zika.tree.geo_inference('region')
        zika.tree.geo_inference('country')
        zika.dump()
        zika.export(prefix='web/data/zika_')
