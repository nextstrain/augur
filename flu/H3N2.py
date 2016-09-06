from __future__ import division, print_function
from collections import defaultdict
from base.process import process
import numpy as np
from datetime import datetime
from base.io_util import myopen

regions = ['SouthAsia', 'Europe', 'China', 'NorthAmerica',
           'China', 'SouthAmerica', 'JapanKorea', 'Oceania' ]


class flu_process(process):
    """process influenza virus sequences in mutliple steps to allow visualization in browser
        * filtering and parsing of sequences
        * alignment
        * tree building
        * frequency estimation of clades and mutations
        * export as json
    """

    def __init__(self, *args, **kwargs):
        super(flu_process, self).__init__(*args, **kwargs)
        print("Initializing flu_process")
        self.HI_titer_fname = self.input_data_path+'_hi_titers.tsv'
        self.HI_strains_fname = self.input_data_path+'_hi_strains.tsv'


    def subsample(self, sampling_threshold, **kwargs):
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

        self.seqs.subsample(category = sampling_category,
                            threshold=sampling_threshold,
                            priority=sampling_priority, **kwargs)
        #tmp = []
        #for seq in self.seqs.seqs.values():
        #    tmp.append((seq.name, sampling_priority(seq), seq.attributes['region'], seq.attributes['date']))
        #print(sorted(tmp, key=lambda x:x[1]))
        #self.seqs.subsample(category = lambda x:(x.attributes['date'].year,x.attributes['date'].month),
        #                    threshold=params.viruses_per_month, repeated=True)




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
        for n in self.tree.tree.find_clades():
            n.attr['cTiter'] = n.cTiter
            n.attr['dTiter'] = n.dTiter

        # SUBSTITUTION MODEL
        self.HI_subs = substitution_model(self.tree.tree, titer_fname = self.HI_titer_fname,**kwargs)
        self.HI_subs.prepare(**kwargs)
        self.HI_subs.train(**kwargs)


    def HI_export(self):
        from base.io_util import write_json
        prefix = self.build_data_path
        if hasattr(self, 'HI_tree'):
            # export the raw titers
            hi_data = self.HI_tree.compile_titers()
            write_json(hi_data, prefix+'titers.json')
            # export the tree model (avidities and potencies only)
            tree_model = {'potency':self.HI_tree.compile_potencies(),
                          'avidity':self.HI_tree.compile_virus_effects(),
                          'dTiter':{n.clade:n.dTiter for n in self.tree.tree.find_clades() if n.dTiter>1e-6}}
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
        node.attr['ep'] = epitope_distance(total_aa_seq, root_total_aa_seq)
        node.attr['ne'] = nonepitope_distance(total_aa_seq, root_total_aa_seq)
        node.attr['rb'] = receptor_binding_distance(total_aa_seq, root_total_aa_seq)


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
    input_data_path = '../nextstrain-db/data/'+lineage
    store_data_path = 'store/'+lineage + '_' + params.resolution +'_'
    build_data_path = 'build/'+lineage + '_' + params.resolution +'_'

    ppy = 12
    flu = flu_process(input_data_path = input_data_path, store_data_path = store_data_path,
                   build_data_path = build_data_path, reference='flu/metadata/h3n2_outgroup.gb',
                   proteins=['SigPep', 'HA1', 'HA2'],
                   method='SLSQP', inertia=np.exp(-1.0/ppy), stiffness=50./ppy)


    flu.load_sequences(fields={0:'strain', 2:'isolate_id', 3:'date', 4:'region',
                         5:'country', 7:"city", 12:"subtype",13:'lineage'})

    time_interval = [datetime.strptime(x, '%Y-%m-%d').date()  for x in params.time_interval]
    pivots = np.arange(time_interval[0].year, time_interval[1].year+time_interval[1].month/12.0, 1.0/ppy)
    flu.seqs.filter(lambda s:
        s.attributes['date']>=time_interval[0] and s.attributes['date']<time_interval[1])

    flu.subsample(params.viruses_per_month)
    flu.align()
    flu.dump()
    # first estimate frequencies globally, then region specific
    flu.estimate_mutation_frequencies(region="global", pivots=pivots)
    for region in regions:
        flu.estimate_mutation_frequencies(region=region)

    flu.subsample(5, repeated=True)
    flu.align()
    flu.build_tree()
    flu.clock_filter(n_iqd=3, plot=True)
    flu.annotate_tree(Tc=0.005, timetree=True, reroot='best')
    flu.tree.geo_inference('region')

    flu.estimate_tree_frequencies()
    flu.dump()

    flu.HI_model()
    H3N2_scores(flu.tree.tree)
    flu.dump()

    flu.export(extra_attr=['serum'])
    flu.HI_export()
