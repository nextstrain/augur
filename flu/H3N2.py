from __future__ import division, print_function
import sys, os, time, gzip, glob
sys.path.append('/ebio/ag-neher/share/users/rneher/nextstrain/nextstrain-augur')
from collections import defaultdict
from base.io_util import make_dir, remove_dir, tree_to_json, write_json, myopen
from base.sequences import sequence_set, num_date
from base.tree import tree
from base.frequencies import alignment_frequencies, tree_frequencies
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import numpy as np
from datetime import datetime

regions = ['SouthAsia', 'Europe', 'China', 'NorthAmerica',
           'China', 'SouthAmerica', 'JapanKorea', 'Oceania' ]

def fix_name(name):
    '''
    function to harmonize names -- will become obsolete once all parsing
    is done by the data base
    '''
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

    def __init__(self, data_path = 'data/h3n2', dt = 0.25, time_interval = ('2008-01-01', '2016-01-01'),
                 out_specs={'out_dir':'data/', 'prefix':'h3n2_', 'qualifier':''},
                 **kwargs):
        super(flu_process, self).__init__()
        print("Initializing flu_process for time interval",time_interval)
        self.data_path = data_path
        self.kwargs = kwargs
        self.out_specs = out_specs
        self.time_interval = [datetime.strptime(time_interval[0], "%Y-%m-%d").date(),
                              datetime.strptime(time_interval[1], "%Y-%m-%d").date()]
        self.filenames()

        # set up containers for frequency estimation and pivots to use
        self.frequencies = defaultdict(dict)
        self.frequency_confidence = defaultdict(dict)
        self.tip_count = defaultdict(dict)
        self.tree_frequencies = defaultdict(dict)
        self.tree_frequency_confidence = defaultdict(dict)
        self.freq_pivot_dt = dt
        self.pivots = np.arange(num_date(self.time_interval[0]),
                                  num_date(self.time_interval[1])+self.freq_pivot_dt,
                                  self.freq_pivot_dt)
        self.seqs=None

        # parse the outgroup information
        if 'outgroup' in kwargs:
            outgroup_file = kwargs['outgroup']
        else:
            outgroup_file = 'flu/metadata/'+out_specs['prefix']+'outgroup.gb'
        self.load_outgroup(outgroup_file)


    def load_outgroup(self,outgroup_file):
        self.outgroup = 'A/Beijing/32/1992' #(tmp_outgroup.features[0].qualifiers['strain'][0]).lower()
        self.outgroup_seq = SeqIO.read(outgroup_file, 'genbank')
        self.genome_annotation = self.outgroup_seq.features

        # grap annotation from genbank
        self.proteins = {f.qualifiers['gene'][0]:FeatureLocation(start=f.location.start, end=f.location.end, strand=1)
                for f in self.genome_annotation if 'gene' in f.qualifiers and f.qualifiers['gene'][0] in ['SigPep', 'HA1', 'HA2']}


    def load_sequences(self):
        # instantiate and population the sequence objects
        self.seqs = sequence_set(self.sequence_fname, reference=self.outgroup)
        self.seqs.ungap()
        self.seqs.parse({0:'strain', 2:'isolate_id', 3:'date', 4:'region',
                         5:'country', 7:"city", 12:"subtype",13:'lineage'}, strip='_')

        # make sure the reference is part of the sequence set
        if self.outgroup in self.seqs.raw_seqs:
            self.seqs.raw_seqs[self.outgroup].seq=self.outgroup_seq.seq
        else:
            print('Outgroup is not in data base')
            self.seqs.raw_seqs[self.outgroup]=self.outgroup_seq

        self.seqs.reference = self.seqs.raw_seqs[self.outgroup]
        # throw out sequences without dates
        self.seqs.parse_date(["%Y-%m-%d"], prune=True)


    def filenames(self):
        '''
        define filenames of input files and intermediates outputs
        '''
        self.HI_strains_fname = self.data_path+'_hi_strains.tsv'
        self.HI_titer_fname = self.data_path+'_hi_titers.tsv'
        self.sequence_fname = self.data_path+'.fasta'

        out_data_path = self.out_specs['out_dir']+self.out_specs['prefix']
        out_data_path += self.out_specs['qualifier']
        self.file_dumps = {}
        self.file_dumps['seqs'] = out_data_path+'sequences.pkl.gz'
        self.file_dumps['tree'] = out_data_path+'tree.newick'
        self.file_dumps['nodes'] = out_data_path+'nodes.pkl.gz'
        self.file_dumps['frequencies'] = out_data_path+'frequencies.pkl.gz'
        self.file_dumps['tree_frequencies'] = out_data_path+'tree_frequencies.pkl.gz'


    def dump(self):
        '''
        write the current state to file
        '''
        from cPickle import dump
        from Bio import Phylo
        for attr_name, fname in self.file_dumps.iteritems():
            if hasattr(self,attr_name):
                print("dumping",attr_name)
                if attr_name=='seqs': self.seqs.raw_seqs = None
                with myopen(fname, 'wb') as ofile:
                    if attr_name=='nodes':
                        continue
                    elif attr_name=='tree':
                        #biopython trees don't pickle well, write as newick + node info
                        self.tree.dump(fname, self.file_dumps['nodes'])
                    else:
                        dump(getattr(self,attr_name), ofile, -1)

    def load(self):
        '''
        reconstruct instance from files
        '''
        from cPickle import load
        for attr_name, fname in self.file_dumps.iteritems():
            if attr_name=='tree':
                continue
            if os.path.isfile(fname):
                with myopen(fname, 'r') as ifile:
                    print('loading',attr_name,'from file',fname)
                    setattr(self, attr_name, load(ifile))

        tree_name = self.file_dumps['tree']
        if os.path.isfile(tree_name):
            if os.path.isfile(self.file_dumps['nodes']):
                node_file = self.file_dumps['nodes']
            else:
                node_file = None
            # load tree, build if no tree file available
            self.build_tree(tree_name, node_file, root='none')


    def subsample(self, sampling_threshold):
        if self.seqs is None:
            self.load_sequences()

        def sampling_category(x):
            return (x.attributes['region'],
                    x.attributes['date'].year,
                    x.attributes['date'].month)

        # load HI titer count to prioritize sequences
        HI_titer_count = {}
        with myopen(self.HI_strains_fname,'r') as ifile:
            for line in ifile:
                strain, count = line.strip().split()
                HI_titer_count[strain]=int(count)

        def sampling_priority(seq):
            sname = seq.attributes['strain']
            if sname in HI_titer_count:
                pr = HI_titer_count[sname]
            else:
                pr = 0
            return pr + len(seq.seq)*0.0001 - 0.01*np.sum([seq.seq.count(nuc) for nuc in 'NRWYMKSHBVD'])

        # constrain to time interval
        self.seqs.raw_seqs = {k:s for k,s in self.seqs.raw_seqs.iteritems() if
                                        s.attributes['date']>=self.time_interval[0] and
                                        s.attributes['date']<self.time_interval[1]}
        self.seqs.subsample(category = sampling_category,
                            threshold=sampling_threshold,
                            priority=sampling_priority )
        #tmp = []
        #for seq in self.seqs.seqs.values():
        #    tmp.append((seq.name, sampling_priority(seq), seq.attributes['region'], seq.attributes['date']))
        #print(sorted(tmp, key=lambda x:x[1]))
        #self.seqs.subsample(category = lambda x:(x.attributes['date'].year,x.attributes['date'].month),
        #                    threshold=params.viruses_per_month, repeated=True)


    def align(self):
        '''
        align sequences, remove non-reference insertions, outlier sequences, and translate
        '''
        self.seqs.align()
        self.seqs.strip_non_reference()
        self.seqs.clock_filter(n_iqd=3, plot=False, max_gaps=0.05, root_seq=self.outgroup)
        self.seqs.translate(proteins=self.proteins)


    def estimate_mutation_frequencies(self, region="global"):
        '''
        calculate the frequencies of mutation in a particular region
        currently the global frequencies should be estimated first
        because this defines the set of positions at which frequencies in
        other regions are estimated.
        '''
        def filter_alignment(aln, region=None, lower_tp=None, upper_tp=None):
            from Bio.Align import MultipleSeqAlignment
            tmp = aln
            if region is not None:
                tmp = [s for s in aln if s.attributes['region']==region]
            if lower_tp is not None:
                tmp = [s for s in aln if s.attributes['num_date']>=lower_tp]
            if upper_tp is not None:
                tmp = [s for s in aln if s.attributes['num_date']<upper_tp]
            return MultipleSeqAlignment(aln)


        if not hasattr(self.seqs, 'aln'):
            print("Align sequences first")
            return

        # loop over nucleotide sequences and translations and calcuate
        # region specific frequencies of mutations above a certain threshold
        for prot, aln in [('nuc',self.seqs.aln)]+ self.seqs.translations.items():
            if region=="global":
                tmp_aln = filter_alignment(aln, lower_tp=self.pivots[0], upper_tp=self.pivots[-1])
                include_set=[]
            else:
                tmp_aln = filter_alignment(aln, region=region, lower_tp=self.pivots[0], upper_tp=self.pivots[-1])
                include_set = set([pos for (pos, mut) in self.frequencies[('global', prot)]])
            time_points = [x.attributes['num_date'] for x in tmp_aln]
            if len(time_points)==0:
                print('no samples in region', region, prot)
                continue

            aln_frequencies = alignment_frequencies(tmp_aln, time_points,
                                            self.pivots, ws=max(2,len(time_points)//10),
                                            **self.kwargs)
            aln_frequencies.mutation_frequencies(min_freq=0.01)
            self.frequencies[(region,prot)] = aln_frequencies.frequencies
            self.frequency_confidence[(region,prot)] = aln_frequencies.calc_confidence()
        self.tip_count[region]=aln_frequencies.counts


    def estimate_tree_frequencies(self, region='global'):
        '''
        estimate frequencies of clades in the tree, possibly region specific
        '''
        if region=='global':
            node_filter_func = None
        else:
            node_filter_func = lambda x:x.region==region

        tree_freqs = tree_frequencies(self.tree.tree, self.pivots,
                                      node_filter = node_filter_func,
                                      ws = max(2,self.tree.tree.count_terminals()//10),
                                      **self.kwargs)

        tree_freqs.estimate_clade_frequencies()
        conf = tree_freqs.calc_confidence()
        self.tree_frequencies[region] = tree_freqs.frequencies
        self.tree_frequency_confidence[region] = conf
        self.tree_pivots = tree_freqs.pivots


    def count_mutations_per_site(self):
        '''
        count the number of independent mutations at each site
        '''
        def mut_struct():
            return defaultdict(int)
        mutation_dict = defaultdict(mut_struct)
        for node in self.tree.tree.find_clades():
            for prot in node.aa_mutations:
                for anc, pos, der in node.aa_mutations[prot]:
                    mutation_dict[(prot,pos)][anc+'->'+der]+=1

        for key in mutation_dict:
            mutation_dict[key]['total'] = np.sum(mutation_dict[key].values())

        self.mutation_count = mutation_dict


    def calculate_associations(self, covariate='passage', lookup=None):
        '''
        calculate the association of amino acid state and
        sequence properties such as passage
        '''
        if not hasattr(self, 'mutation_count'):
            self.count_mutations_per_site()

        # calculate associations
        from scipy.stats import chi2_contingency
        self.associations = {}
        if lookup is None:
            lookup=lambda x:x

        # loop over all positions (currently rather clumsy)
        for prot, pos in self.mutation_count:
            assoc = defaultdict(int)
            for node in self.tree.tree.get_terminals(): # extract info from each node
                if hasattr(node, covariate):
                    assoc[(node.translations[prot][pos-1], lookup(node.__getattribute__(covariate)))]+=1

            # make contingency matrix
            aa_states = sorted(set([x[0] for x in assoc]))
            cov_states = sorted(set([x[1] for x in assoc]))
            contingeny_matrix = np.zeros((len(aa_states), len(cov_states)))
            for a, c in assoc:
                contingeny_matrix[aa_states.index(a), cov_states.index(c)] = assoc[(a,c)]
            print(contingeny_matrix, assoc)
            g, p, dof, expctd = chi2_contingency(contingeny_matrix, lambda_="log-likelihood")
            assoc['contingency matrix'] = contingeny_matrix
            assoc['aa']=aa_states
            assoc['covariates']=cov_states
            assoc['g_test'] = (g,p)

            self.associations[(prot, pos)] = assoc


    def build_tree(self, infile=None, nodefile=None, root='best', Tc=0.01):
        '''
        instantiate a tree object and make a time tree
        if infiles are None, the tree is build from scratch. Otherwise
        the tree is loaded from file
        '''
        self.tree = tree(aln=self.seqs.aln, proteins = self.proteins)
        if infile is None:
            self.tree.build(root=root)
        else:
            self.tree.tt_from_file(infile, nodefile=nodefile, root=root)
        # if node file is none, no time information is available.
        # hence make a coalescent model time tree and decorate that tree.
        if self.rm_og:
            flu.remove_outgroup()
            flu.tree.tt.prepare_tree()

        if nodefile is None:
            self.tree.timetree(Tc=Tc, infer_gtr=True)
            self.tree.add_translations()
            self.tree.refine()
            self.tree.layout()


    def remove_outgroup(self):
        '''
        THIS ASSUMES LADDERIZATION!
        '''
        while np.sum([x.branch_length for x in self.tree.tree.root.clades])>0.02:
            self.tree.tree.root = self.tree.tree.root.clades[1]
            self.tree.tree.root.branch_length=0.001
            self.tree.tree.root.up = None


    def HI_model(self, **kwargs):
        '''
        estimate a tree and substitution model using titers titer_fname.
        '''
        from base.titer_model import tree_model, substitution_model
        ## TREE MODEL
        self.HI_tree = tree_model(self.tree.tree, titer_fname = self.HI_titer_fname, **kwargs)
        self.HI_tree.prepare(**kwargs)
        self.HI_tree.train(**kwargs)
        # add tree attributes to the list of attributes that are saved in intermediate files
        self.tree.dump_attr.extend(['cTiter', 'dTiter'])

        # SUBSTITUTION MODEL
        self.HI_subs = substitution_model(self.tree.tree, titer_fname = self.HI_titer_fname,**kwargs)
        self.HI_subs.prepare(**kwargs)
        self.HI_subs.train(**kwargs)


    def export(self, prefix='web/data/', extra_attr = []):
        '''
        export the tree, sequences, frequencies to json files for visualization
        in the browser
        '''

        # export json file that contains alignment diversity column by column
        self.seqs.export_diversity(prefix+'entropy.json')
        # exports the tree and the sequences inferred for all clades in the tree
        self.tree.export(path=prefix, extra_attr = extra_attr + ["subtype", "country", "region", "nuc_muts",
                            "ep", "ne", "rb", "aa_muts","lab", "accession","isolate_id"])


        # local function or round frequency estimates to useful precision (reduces file size)
        def process_freqs(freq):
            return [round(x,4) for x in freq]

        # construct a json file containing all frequency estimate
        # the format is region_protein:159F for mutations and region_clade:123 for clades
        freq_json = {'pivots':process_freqs(self.pivots)}
        freq_json['counts'] = {x:list(counts) for x, counts in self.tip_count.iteritems()}
        for (region, gene), tmp_freqs in self.frequencies.iteritems():
            for mut, freq in tmp_freqs.iteritems():
                label_str =  region+"_"+ gene + ':' + str(mut[0]+1)+mut[1]
                freq_json[label_str] = process_freqs(freq)
        # repeat for clade frequencies in trees
        for region in self.tree_frequencies:
            for clade, freq in self.tree_frequencies[region].iteritems():
                label_str = region+'_clade:'+str(clade)
                freq_json[label_str] = process_freqs(freq)
        # write to one frequency json
        write_json(freq_json, prefix+'frequencies.json', indent=None)

        self.HI_export(prefix)

    def HI_export(self, prefix):
        if hasattr(self, 'HI_tree'):
            # export the raw titers
            hi_data = self.HI_tree.compile_titers()
            write_json(hi_data, prefix+'titers.json')
            # export the tree model (avidities and potencies only)
            tree_model = {'potency':self.HI_tree.compile_potencies(),
                          'avidities':self.HI_tree.compile_virus_effects()}
            write_json(tree_model, prefix+'titer_tree_model.json')
        else:
            print('Tree model not yet trained')

        if hasattr(self, 'HI_tree'):
            # export the substitution model
            subs_model = {'potency':self.HI_subs.compile_potencies(),
                          'avidity':self.HI_subs.compile_virus_effects(),
                          'substitution':self.HI_subs.compile_substitution_effects()}
            write_json(subs_model, prefix+'titer_subs_model.json')
        else:
            print('Substitution model not yet trained')


def H3N2_scores(tree, epitope_mask_version='wolf'):
    '''
    takes a H3N2 HA tree and assigns H3 specific characteristics to
    internal and external nodes
    '''
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
    with open('flu/metadata/h3n2_epitope_masks.tsv') as f:
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


def plot_frequencies(flu, plot_regions):
    pivots = flu.pivots
    plt.figure()


def plot_trace(ax, pivots, freq, err, n_std_dev=1, err_smoothing=3, show_errorbars=True, c='r', ls='-', label=None):
    ax.plot(pivots, freq, c=c, ls=ls, label=label)
    if show_errorbars:
        smerr = 1.0/np.convolve(1.0/err, np.ones(err_smoothing, dtype=float)/err_smoothing, mode='same')
        ax.fill_between(pivots, np.maximum(0,freq-n_std_dev*smerr),
                np.minimum(1,freq+n_std_dev*smerr),
                facecolor=c, linewidth=0, alpha=0.1)

def plot_frequencies(flu, gene, mutation=None, plot_regions=None, all_muts=False, **kwargs):
    import seaborn as sns
    sns.set_style('whitegrid')
    cols = sns.color_palette()
    linestyles = ['-', '--', '-.', ':']
    if plot_regions is None:
        plot_regions=regions
    pivots = flu.pivots
    plt.figure()
    ax=plt.subplot(111)
    if type(mutation)==int:
        mutations = [x for x,freq in flu.frequencies[(gene, 'global')].iteritems()
                     if (x[0]==mutation)&(freq[0]<0.5 or all_muts)]
    elif mutation is not None:
        mutations = [mutation]
    else:
        mutations=None

    if mutations is None:
        for ri, region in enumerate(plot_regions):
            count=flu.tip_count[region]
            plt.plot(pivots, count, c=cols[ri%len(cols)], label=region)
    else:
        print("plotting mutations", mutations)
        for ri,region in enumerate(plot_regions):
            for mi,mut in enumerate(mutations):
                if mut in flu.frequencies[(gene, region)]:
                    freq = flu.frequencies[(gene, region)][mut]
                    err = flu.frequency_confidence[(gene, region)][mut]
                    c=cols[ri%len(cols)]
                    label_str = str(mut[0]+1)+mut[1]+', '+region
                    plot_trace(ax, pivots, freq, err, c=c,
                        ls=linestyles[mi%len(linestyles)],label=label_str, **kwargs)
                else:
                    print(mut, 'not found in region',region)
    ax.ticklabel_format(useOffset=False)
    ax.legend(loc=2)


if __name__=="__main__":
    import argparse
    import matplotlib.pyplot as plt
    plt.ion()

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-v', '--viruses_per_month', type = int, default = 10, help='number of viruses sampled per month')
    parser.add_argument('-y', '--resolution', type = str, default = '3y', help='outfile suffix')
    parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
    parser.add_argument('-d', '--download', action='store_true', default = False, help='load from database')
    parser.add_argument('-t', '--time_interval', nargs=2, default=('2012-01-01', '2016-01-01'),
                            help='time interval to sample sequences from: provide dates as YYYY-MM-DD')
    parser.add_argument('--load', action='store_true', help = 'recover from file')
    params = parser.parse_args()
    lineage = 'h3n2'
    data_path = '../nextstrain-db/data/'+lineage

    out_specs = {'out_dir':'data/',
                 'prefix':lineage+'_', 'qualifier':params.resolution+'_'}

    ppy = 4
    flu = flu_process(method='SLSQP', dtps=2.0, stiffness=50./ppy, dt=1.0/ppy,
                      time_interval=params.time_interval,
                      inertia=np.exp(-1.0/ppy), data_path = data_path, out_specs=out_specs)
    remove_outgroup = int(params.time_interval[0][:4])>1996
    flu.rm_og = remove_outgroup
    if params.load:
        flu.load()
        H3N2_scores(flu.tree.tree)
    else:
        remove_outgroup = int(params.time_interval[0][:4])>1996
        flu.load_sequences()
        flu.subsample(params.viruses_per_month)
        flu.align()
        flu.dump()
        # first estimate frequencies globally, then region specific
        flu.estimate_mutation_frequencies(region="global")
        #for region in regions:
        #    flu.estimate_mutation_frequencies(region=region)
        flu.dump()
        flu.build_tree(Tc=0.005)

        flu.tree.add_translations()
        flu.tree.refine()
        flu.tree.layout()
        flu.tree.geo_inference('region')
        #flu.tree.geo_inference('country')
        flu.dump()
        flu.estimate_tree_frequencies()
        flu.dump()

        flu.HI_model()
        H3N2_scores(flu.tree.tree)
        flu.dump()

        flu.export(prefix='json/'+out_specs['prefix']+out_specs['qualifier'],
                   extra_attr=['cTiter', 'dTiter', 'aa_mut_str', 'serum'])
