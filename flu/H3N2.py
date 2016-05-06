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
    fields = tmp_name.split('/')
    # fix two digit dates
    if len(fields[-1])==2:
        try:
            y = int(fields[-1])
            if y>16:
                y=1900+y
            else:
                y=2000+y
            new_name =  '/'.join(fields[:-1])+'/'+str(y)
        except:
            new_name = tmp_name
    else:
        new_name = tmp_name
    return new_name.strip().strip('_')

class flu_process(object):
    """process influenza virus sequences in mutliple steps to allow visualization in browser
        * filtering and parsing of sequences
        * alignment
        * tree building
        * frequency estimation of clades and mutations
        * export as json
    """

    def __init__(self, fname = 'data/H3N2_gisaid_epiflu_sequence.fasta',
                 out_specs={'data_dir':'data/', 'prefix':'H3N2_', 'qualifier':''},
                 **kwargs):
        super(flu_process, self).__init__()
        self.fname = fname
        self.kwargs = kwargs
        self.out_specs = out_specs
        if 'outgroup' in kwargs:
            outgroup_file = kwargs['outgroup']
        else:
            outgroup_file = 'source_data/'+out_specs['prefix']+'outgroup.gb'
        tmp_outgroup = SeqIO.read(outgroup_file, 'genbank')
        self.outgroup = tmp_outgroup.features[0].qualifiers['strain'][0]
        genome_annotation = tmp_outgroup.features
        ref_seq = SeqIO.read(outgroup_file, 'genbank')
        self.proteins = {f.qualifiers['gene'][0]:FeatureLocation(start=f.location.start, end=f.location.end, strand=1)
                for f in ref_seq.features if 'gene' in f.qualifiers and f.qualifiers['gene'][0] in ['SigPep', 'HA1', 'HA2']}

        self.time_interval = [datetime.strptime('2008-01-01', "%Y-%m-%d").date(),
                              datetime.strptime('2016-01-01', "%Y-%m-%d").date()]
        self.frequencies = defaultdict(dict)
        self.pivots = np.linspace(num_date(self.time_interval[0]),
                                  num_date(self.time_interval[1]),40)

        self.seqs = sequence_set(self.fname, reference=self.outgroup)
        self.seqs.ungap()
        self.seqs.parse({0:'strain', 1:'isolate_id', 3:'passage', 5:'date', 7:'lab', 8:"accession"}, strip='_')
        self.fix_strain_names()
        self.seqs.raw_seqs[self.outgroup].seq=tmp_outgroup.seq
        self.seqs.raw_seqs['A/Beijing/32/1992'].attributes['date']='1992-01-01'
        self.seqs.reference = self.seqs.raw_seqs[self.outgroup]
        self.seqs.parse_date(["%Y-%m-%d"], prune=True)
        self.geo_parse()
        self.filenames()


    def filenames(self):
        data_path = self.out_specs['data_dir']+self.out_specs['prefix']
        self.HI_strains_fname = data_path+'HI_strains.txt'
        self.HI_titer_fname = data_path+'HI_titers.txt'
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


    def geo_parse(self):
        import csv,re
        reader = csv.DictReader(open("source_data/geo_synonyms.tsv"), delimiter='\t')       # list of dicts
        label_to_country = {}
        for line in reader:
            label_to_country[line['label'].lower()] = line['country']
        for strain, v in self.seqs.raw_seqs.iteritems():
            if "country" not in v.attributes:
                v.attributes['country'] = 'Unknown'
                try:
                    fixed_strain = strain.replace('_',' ')
                    label = re.match(r'^[AB]/([^/]+)/', fixed_strain).group(1).lower()                       # check first for whole geo match
                    if label in label_to_country:
                        v.attributes['country'] = label_to_country[label]
                    else:
                        label = re.match(r'^[AB]/([^\-^\/]+)[\-\/]', fixed_strain).group(1).lower()          # check for partial geo match
                    if label in label_to_country:
                        v.attributes['country'] = label_to_country[label]
                    else:
                        label = re.match(r'^[AB]/([A-Z][a-z]+)[A-Z0-9]', fixed_strain).group(1).lower()          # check for partial geo match
                    if label in label_to_country:
                        v.attributes['country'] = label_to_country[label]
                    if v.attributes['country'] == 'Unknown':
                        print("couldn't parse location for", fixed_strain)
                except:
                    print("couldn't parse location for", fixed_strain, label)

        reader = csv.DictReader(open("source_data/geo_regions.tsv"), delimiter='\t')        # list of dicts
        country_to_region = {}
        for line in reader:
            country_to_region[line['country']] = line['region']
        for strain, v in self.seqs.raw_seqs.iteritems():
            v.attributes['region'] = 'Unknown'
            if v.attributes['country'] in country_to_region:
                v.attributes['region'] = country_to_region[v.attributes['country']]
            if v.attributes['country'] != 'Unknown' and v.attributes['region'] == 'Unknown':
                print("couldn't parse region for", strain, "country:", v.attributes["country"])


    def subsample(self):
        HI_titer_count = {}
        with myopen(self.HI_strains_fname,'r') as ifile:
            for line in ifile:
                strain, count = line.strip().split()
                HI_titer_count[strain]=int(count)

        def sampling_priority(seq):
            sname = seq.attributes['strain'].upper()
            if sname in HI_titer_count:
                pr = HI_titer_count[sname]
            else:
                pr = 0
            return pr + len(seq.seq)*0.0001 - 0.01*np.sum([seq.seq.count(nuc) for nuc in 'NRWYMKSHBVD'])

        self.seqs.raw_seqs = {k:s for k,s in self.seqs.raw_seqs.iteritems() if
                                        s.attributes['date']>=self.time_interval[0] and
                                        s.attributes['date']<self.time_interval[1]}
        self.seqs.subsample(category = lambda x:(x.attributes['region'],
                                                 x.attributes['date'].year,
                                                 x.attributes['date'].month),
                            threshold=params.viruses_per_month, priority=sampling_priority )
        #tmp = []
        #for seq in self.seqs.seqs.values():
        #    tmp.append((seq.name, sampling_priority(seq), seq.attributes['region'], seq.attributes['date']))
        #print(sorted(tmp, key=lambda x:x[1]))
        #self.seqs.subsample(category = lambda x:(x.attributes['date'].year,x.attributes['date'].month),
        #                    threshold=params.viruses_per_month, repeated=True)


    def align(self):
        self.seqs.align()
        self.seqs.strip_non_reference()
        self.seqs.clock_filter(n_iqd=3, plot=True, max_gaps=0.05, root_seq=self.outgroup)
        self.seqs.translate(proteins=self.proteins)

    def estimate_mutation_frequencies(self):
        if not hasattr(self.seqs, 'aln'):
            print("Align sequences first")
            return
        time_points = [x.attributes['num_date'] for x in self.seqs.aln]
        aln_frequencies = alignment_frequencies(self.seqs.aln, time_points,
                                self.pivots, ws=len(time_points)/10, **self.kwargs)
        aln_frequencies.mutation_frequencies(min_freq=0.1)
        self.frequencies['nuc'] = aln_frequencies.frequencies
        for prot in self.seqs.translations:
            aln_frequencies = alignment_frequencies(self.seqs.translations[prot], time_points,
                                            self.pivots, ws=len(time_points)//10, **self.kwargs)
            aln_frequencies.mutation_frequencies(min_freq=0.01)
            self.frequencies[prot] = aln_frequencies.frequencies

    def estimate_tree_frequencies(self):
        tree_freqs = tree_frequencies(self.tree.tree, self.pivots,
                                      ws = self.tree.tree.count_terminals()//10, **self.kwargs)
        tree_freqs.estimate_clade_frequencies()
        self.tree_frequencies = tree_freqs.frequencies


    def build_tree(self, infile=None):
        self.tree = tree(aln=self.seqs.aln, proteins = self.proteins)
        if infile is None:
            self.tree.build()
        else:
            self.tree.tt_from_file(infile)
        self.tree.timetree(Tc=0.01, infer_gtr=True)
        self.tree.add_translations()
        self.tree.refine()
        self.tree.layout()


    def export(self, prefix='web/data/'):
        def process_freqs(freq):
            return [round(x,4) for x in freq]
        self.seqs.export_diversity(prefix+'entropy.json')
        self.tree.export(path=prefix, extra_attr = ["subtype", "country", "region", "nuc_muts",
                                                    "ep", "ne", "rb", "aa_muts","lab", "accession","isolate_id"])

        freq_json = {'pivots':process_freqs(self.pivots)}
        for gene, tmp_freqs in self.frequencies.iteritems():
            for mut, freq in tmp_freqs.iteritems():
                freq_json['_'.join([gene, str(mut[0]+1), mut[1]])] = process_freqs(freq)
        for clade, freq in self.tree_frequencies.iteritems():
            freq_json['clade_'+str(clade)] = process_freqs(freq)
        write_json(freq_json, prefix+'frequencies.json', indent=None)


    def HI_model(self):
        from HI_model import HI_tree
        self.HI = HI_tree(self.tree.tree, HI_fname = self.HI_titer_fname)
        self.HI.map_HI(training_fraction=1.0,method='nnl1reg', map_to_tree=True)
        self.HI.add_titers()


def H3N2_scores(tree, epitope_mask_version='wolf'):
    def epitope_sites(aa):
        return aa[epitope_mask[:len(aa)]]

    def nonepitope_sites(aa):
        return aa[~epitope_mask[:len(aa)]]

    def receptor_binding_sites(aa):
        '''
        Receptor binding site mutations from Koel et al. 2014
        These are (145, 155, 156, 158, 159, 189, 193) in canonical HA numbering
        need to subtract one since python arrays start at 0
        '''
        sp = 16
        rbs = map(lambda x:x+sp-1, [145, 155, 156, 158, 159, 189, 193])
        return np.array([aa[pos] for pos in rbs])

    def get_total_peptide(node):
        '''
        the concatenation of signal peptide, HA1, HA1
        '''
        return np.fromstring(node.translations['SigPep']+node.translations['HA1']
                           + node.translations['HA2'], 'S1')

    def epitope_distance(aaA, aaB):
        """Return distance of sequences aaA and aaB by comparing epitope sites"""
        epA = epitope_sites(aaA)
        epB = epitope_sites(aaB)
        distance = np.sum(epA!=epB)
        return distance

    def nonepitope_distance(aaA, aaB):
        """Return distance of sequences aaA and aaB by comparing non-epitope sites"""
        neA = nonepitope_sites(aaA)
        neB = nonepitope_sites(aaB)
        distance = np.sum(neA!=neB)
        return distance

    def receptor_binding_distance(aaA, aaB):
        """Return distance of sequences aaA and aaB by comparing receptor binding sites"""
        neA = receptor_binding_sites(aaA)
        neB = receptor_binding_sites(aaB)
        distance = np.sum(neA!=neB)
        return distance

    epitope_map = {}
    with open('source_data/H3N2_epitope_masks.tsv') as f:
        for line in f:
            (key, value) = line.strip().split()
            epitope_map[key] = value
    if epitope_mask_version in epitope_map:
        epitope_mask = np.fromstring(epitope_map[epitope_mask_version], 'S1')=='1'
    root = tree.root
    root_total_aa_seq = get_total_peptide(root)
    for node in tree.find_clades():
        total_aa_seq = get_total_peptide(node)
        node.ep = epitope_distance(total_aa_seq, root_total_aa_seq)
        node.ne = nonepitope_distance(total_aa_seq, root_total_aa_seq)
        node.rb = receptor_binding_distance(total_aa_seq, root_total_aa_seq)



if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-v', '--viruses_per_month', type = int, default = 10, help='number of viruses sampled per month')
    parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
    parser.add_argument('--load', action='store_true', help = 'recover from file')

    params = parser.parse_args()

    flu = flu_process(method='SLSQP', dtps=2.0, stiffness=20, inertia=0.9)
    if params.load:
        flu.load()
    else:
        flu.subsample()
        flu.align()
        flu.dump()
        flu.estimate_mutation_frequencies()
        flu.dump()
        flu.build_tree()
        flu.dump()
        flu.estimate_tree_frequencies()

        H3N2_scores(flu.tree.tree)
        flu.export()
