from __future__ import division, print_function
from collections import defaultdict
import sys
sys.path.append('')  # need to import from base
from base.process import process
import numpy as np
from datetime import datetime, timedelta
from base.io_util import myopen
from Bio import SeqIO
from glob import glob
import os

regions = ['africa', 'south_asia', 'europe', 'china', 'north_america',
           'china', 'south_america', 'japan_korea', 'oceania', 'southeast_asia']

region_groups = {'NA':'north_america',
                 'AS':['china', 'japan_korea', 'south_asia', 'southeast_asia'],
                 'OC':'oceania', 'EU':'europe'}

attribute_nesting = {'geographic location':['region', 'country', 'city'],}

color_options = {
    "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
    "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
    "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
    "ep":{"key":"ep", "legendTitle":"Epitope Mutations", "menuItem":"epitope mutations", "type":"continuous"},
    "ne":{"key":"ne", "legendTitle":"Non-epitope Mutations", "menuItem":"nonepitope mutations", "type":"continuous"},
    "rb":{"key":"rb", "legendTitle":"Receptor Binding Mutations", "menuItem":"RBS mutations", "type":"continuous"},
    "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"}}
panels = ['tree', 'entropy', 'frequencies']


class dengue_process(process):

    def __init__(self, *args, **kwargs):
        super(dengue_process, self).__init__(*args, **kwargs)
        print("Initializing dengue_process")
        self.titer_fname = kwargs['titer_fname']
        self.strains_fname = kwargs['strain_fname']


    def subsample(self, sampling_threshold, **kwargs):

        def sampling_category(x):
            '''
            Subsample per region, per month.
            '''
            return (x.attributes['region'],
                    x.attributes['date'].year,
                    x.attributes['date'].month)

        # load titer count to prioritize sequences
        titer_count = {}
        forced_strains = []
        with myopen(self.strains_fname,'r') as ifile:
            for line in ifile:
                strain, count = line.strip().split()
                titer_count[strain]=int(count)
                if strain.startswith('DENV%s'%serotype):
                    forced_strains.append(strain)

        def sampling_priority(seq):
            '''
            Prefers more titer measurements, longer sequences. Penalizes ambiguity codes.
            '''
            sname = seq.attributes['strain']
            if sname in titer_count:
                pr = titer_count[sname]
            else:
                pr = 0
            return pr + len(seq.seq)*0.0001 - 0.01*np.sum([seq.seq.count(nuc) for nuc in 'NRWYMKSHBVD'])

        self.seqs.subsample(category = sampling_category,
                            threshold=sampling_threshold,
                            priority=sampling_priority,
                            forced_strains=forced_strains,
                            **kwargs)

    def titer_model(self, **kwargs):
        '''
        estimate a titer tree and substitution model using titers in titer_fname.
        '''
        from base.titer_model import tree_model, substitution_model
        ## TREE MODEL
        self.titer_tree = tree_model(self.tree.tree, titer_fname = self.titer_fname, **kwargs)
        self.titer_tree.prepare(**kwargs) # make training set, find subtree with titer measurements, and make_treegraph
        self.titer_tree.train(**kwargs) # pick longest branch on path between each (test, ref) pair, assign titer drops to this branch
                                 # then calculate a cumulative antigenic evolution score for each node
        # add tree attributes to the list of attributes that are saved in intermediate files
        for n in self.tree.tree.find_clades():
            n.attr['cTiter'] = n.cTiter
            n.attr['dTiter'] = n.dTiter

        # SUBSTITUTION MODEL
        self.titer_subs = substitution_model(self.tree.tree, titer_fname = self.titer_fname,**kwargs)
        self.titer_subs.prepare(**kwargs)
        self.titer_subs.train(**kwargs)


    def titer_export(self):
        from base.io_util import write_json
        prefix = self.build_data_path
        if hasattr(self, 'titer_tree'):
            # export the raw titers
            data = self.titer_tree.compile_titers()
            write_json(data, prefix+'titers.json', indent=1)
            # export the tree model (avidities and potencies only)
            tree_model = {'potency':self.titer_tree.compile_potencies(),
                          'avidity':self.titer_tree.compile_virus_effects(),
                          'dTiter':{n.clade:n.dTiter for n in self.tree.tree.find_clades() if n.dTiter>1e-6}}
            write_json(tree_model, prefix+'tree_model.json')
        else:
            print('Tree model not yet trained')

        if hasattr(self, 'titer_tree'):
            # export the substitution model
            titer_subs_model = {'potency':self.titer_subs.compile_potencies(),
                          'avidity':self.titer_subs.compile_virus_effects(),
                          'substitution':self.titer_subs.compile_substitution_effects()}
            write_json(titer_subs_model, prefix+'titer_subs_model.json')
        else:
            print('Substitution model not yet trained')

def plot_trace(ax, pivots, freq, err, n_std_dev=1, err_smoothing=3, show_errorbars=True, c='r', ls='-', label=None):
    ax.plot(pivots, freq, c=c, ls=ls, label=label)
    if show_errorbars:
        smerr = 1.0/np.convolve(1.0/err, np.ones(err_smoothing, dtype=float)/err_smoothing, mode='same')
        ax.fill_between(pivots, np.maximum(0,freq-n_std_dev*smerr),
                np.minimum(1,freq+n_std_dev*smerr),
                facecolor=c, linewidth=0, alpha=0.1)

def plot_frequencies(dengue, gene, mutation=None, plot_regions=None, all_muts=False, ax=None, **kwargs):
    import seaborn as sns
    sns.set_style('whitegrid')
    cols = sns.color_palette()
    linestyles = ['-', '--', '-.', ':']
    if plot_regions is None:
        plot_regions=regions
    pivots = dengue.pivots
    if ax is None:
        plt.figure()
        ax=plt.subplot(111)
    if type(mutation)==int:
        mutations = [x for x,freq in dengue.mutation_frequencies[('global', gene)].iteritems()
                     if (x[0]==mutation)&(freq[0]<0.5 or all_muts)]
    elif mutation is not None:
        mutations = [mutation]
    else:
        mutations=None

    if mutations is None:
        for ri, region in enumerate(plot_regions):
            count=dengue.mutation_frequency_counts[region]
            plt.plot(pivots, count, c=cols[ri%len(cols)], label=region)
    else:
        print("plotting mutations", mutations)
        for ri,region in enumerate(plot_regions):
            for mi,mut in enumerate(mutations):
                if mut in dengue.mutation_frequencies[(region, gene)]:
                    freq = dengue.mutation_frequencies[(region, gene)][mut]
                    err = dengue.mutation_frequency_confidence[(region, gene)][mut]
                    c=cols[ri%len(cols)]
                    label_str = str(mut[0]+1)+mut[1]+', '+region
                    plot_trace(ax, pivots, freq, err, c=c,
                        ls=linestyles[mi%len(linestyles)],label=label_str, **kwargs)
                else:
                    print(mut, 'not found in region',region)
    ax.ticklabel_format(useOffset=False)
    ax.legend(loc=2)

def select_serotype(infile, path, serotype):
    '''
    Takes any-serotype fasta file, path to save output, and desired serotype.
    Writes appropriate subset of sequences to dengue_serotype.fasta
    Returns path to output file as string.
    '''
    firstfield = 'DENV%s'%serotype
    sequences = [ i for i in SeqIO.parse(infile, 'fasta') if i.description.split('/')[0] == firstfield ]
    SeqIO.write(sequences, path+'dengue_%s.fasta'%serotype, 'fasta')
    return path+'dengue_%s.fasta'%serotype

if __name__=="__main__":
    import argparse
    import matplotlib.pyplot as plt
    plt.ion()

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-v', '--viruses_per_month', type = int, default = 10, help='number of viruses sampled per month')
    parser.add_argument('-y', '--years_back', type = str, default = 3, help='number of years back to sample')
    parser.add_argument('--resolution', type = str, help ="outfile suffix, can determine -v and -y")
    parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
    parser.add_argument('-d', '--download', action='store_true', default = False, help='load from database')
    parser.add_argument('-t', '--time_interval', nargs=2, help='specify time interval rather than use --years_back')
    parser.add_argument('-s', '--serotype', type = str, choices=['1', '2', '3', '4', 'any', 'all'], default='any', help = 'which serotype of dengue to build trees for; any = include all serotypes in one build; all = run all 5 builds')
    parser.add_argument('--load', action='store_true', help = 'recover from file')
    parser.add_argument('--no_tree', default=False, action='store_true', help = "don't build a tree")
    parser.add_argument('-tf', '--titer_fname', type = str, default='dengue', help = 'prefix of titer and strain fnames found in fauna/data. E.g., `agm_dengue` --> `agm_dengue_titers.tsv`, `agm_dengue_strains.tsv`')
    params = parser.parse_args()

    # default values for --viruses_per_month and --years_back from resolution
    if params.resolution == "2y":
		params.viruses_per_month = 15
		params.years_back = 2
    elif params.resolution == "3y":
		params.viruses_per_month = 7
		params.years_back = 3
    elif params.resolution == "6y":
		params.viruses_per_month = 3
		params.years_back = 6
    elif params.resolution == "12y":
		params.viruses_per_month = 2
		params.years_back = 12

    # construct time_interval from years_back
    if not params.time_interval:
        today_str = "{:%Y-%m-%d}".format(datetime.today())
        date_str = "{:%Y-%m-%d}".format(datetime.today() - timedelta(days=365.25 * params.years_back))
        params.time_interval = [date_str, today_str]

    # For any-serotype build, use dengue 3 outgroup and look for files like dengue.fasta
    serotype = params.serotype

    if serotype == 'any':
        lineage = 'dengue'
        reference_fname = './dengue/metadata/dengue_%s_outgroup.gb'%'3'
        newest_sequence_file = sorted(glob('../fauna/data/%s.fasta'%lineage), key=lambda f: os.path.getmtime(f))[-1]

    # For serotype-specific build, use the corresponding outgroup
    else:
        lineage = 'dengue_%s'%serotype
        reference_fname = './dengue/metadata/dengue_%s_outgroup.gb'%serotype

        try: # Look for a serotype-specific fasta
            newest_sequence_file = sorted(glob('../fauna/data/%s*.fasta'%lineage), key=lambda f: os.path.getmtime(f))[-1]
        except: # If it doesn't exist, try to pull serotype-specific sequences out of the any-serotype fasta (warn the user of this behavior)
            newest_sequence_file = select_serotype('../fauna/data/dengue.fasta', '../fauna/data/', serotype)
            print('WARNING: Did not find serotype-specific fasta file.\nPulled sequences with serotype %s from any-serotype fasta file %s\nWrote these to file %s'%(serotype, '../fauna/data/dengue.fasta', newest_sequence_file))

    input_data_path = newest_sequence_file.split('.fasta')[0]
    sequence_fname = newest_sequence_file
    titer_fname = '../fauna/data/'+params.titer_fname+'_titers.tsv'
    strain_fname = '../fauna/data/'+params.titer_fname+'_strains.tsv'

    if params.resolution:
        store_data_path = 'store/'+lineage + '_' + params.resolution +'_'
        build_data_path = 'build/'+lineage + '_' + params.resolution +'_'
    else:
        store_data_path = 'store/'+lineage + '_'
        build_data_path = 'build/'+lineage + '_'

    ppy = 12

    dengue = dengue_process(input_data_path = input_data_path, store_data_path = store_data_path,
                   build_data_path = build_data_path, reference='dengue/metadata/'+lineage+'_outgroup.gb',
                   proteins= ['C', 'M', 'E', 'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', '2K', 'NS4B', 'NS5'],
                   method='SLSQP', inertia=np.exp(-1.0/ppy), stiffness=2.*ppy, titer_fname = titer_fname, strain_fname = strain_fname)


    if params.load:
        dengue.load()
        # dengue.export()
    else:
        dengue.load_sequences(fields={0:'strain', 1:'accession', 2:'date', 3:'region', 4:'country',
                        5:'division', 6: 'location'})

        time_interval = [datetime.strptime(x, '%Y-%m-%d').date() for x in params.time_interval]
        # pivots = np.arange(time_interval[0].year+(time_interval[0].month-1)/12.0,
        #                    time_interval[1].year+time_interval[1].month/12.0, 1.0/ppy)

        dropped_strains = ['DENV1/VIETNAM/BIDV992/2006', 'DENV1/FRANCE/00475/2008', 'DENV1/VIETNAM/BIDV3990/2008', 'DENV2/HAITI/DENGUEVIRUS2HOMOSAPIENS1/2016',# probable recombinants
        'DENV2/AUSTRALIA/QML22/2015' # Wrong serotype? Abnormal date.
        ]
        dengue.seqs.filter(lambda s: s.id not in dropped_strains)
        dengue.seqs.filter(lambda s:
            s.attributes['date']>=time_interval[0] and s.attributes['date']<time_interval[1])
        print(len(dengue.seqs.all_seqs))
        dengue.seqs.filter(lambda s: len(s.seq)>=900)
        dengue.subsample(params.viruses_per_month)
        dengue.align(debug=True)

        dengue.dump()
        # first estimate frequencies globally, then region specific
        # dengue.estimate_mutation_frequencies(region="global", pivots=pivots)
        # for region in region_groups.iteritems():
        #     dengue.estimate_mutation_frequencies(region=region)

        if not params.no_tree:
            dengue.align(debug=True)
            dengue.build_tree(debug=True)
            dengue.clock_filter(n_iqd=3, plot=True)
            dengue.annotate_tree(Tc=0.005, timetree=True, reroot='best')
            # dengue.tree.geo_inference('region')
            # dengue.estimate_tree_frequencies()
            dengue.dump()

    dengue.titer_model()
    # H3N2_scores(dengue.tree.tree)
    dengue.dump()
    # dengue.matchClades(clade_designations[lineage])
    dengue.titer_export()
    dengue.export(extra_attr=['serum'], controls=attribute_nesting, color_options=color_options, panels=panels )
