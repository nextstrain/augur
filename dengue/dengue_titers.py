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

            if seq.id in forced_strains:
                pr += 10.0

            return pr + len(seq.seq)*0.00005 - 0.01*np.sum([seq.seq.count(nuc) for nuc in 'NRWYMKSHBVD'])

        self.seqs.subsample(category = sampling_category,
                            threshold=sampling_threshold,
                            priority=sampling_priority,
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

        if kwargs['training_fraction'] != 1.0:
            self.titer_tree.validate() #(plot=True, fname='treeModel_%s.png'%lineage)
            self.titer_subs.validate() #(plot=True, fname='subModel_%s.png'%lineage)
            logfile = open('params_log.tsv', 'a')
            logfile.write('\t'.join(['substitution_model', lineage, '%.3f'%kwargs['lam_pot'], '%.3f'%kwargs['lam_avi'], '%.3f'%kwargs['lam_drop'], '%.3f'%self.titer_subs.abs_error, '%.3f'%self.titer_subs.rms_error, '%.3f'%self.titer_subs.r2])+'\n')
            logfile.write('\t'.join(['tree_model', lineage, '%.3f'%kwargs['lam_pot'], '%.3f'%kwargs['lam_avi'], '%.3f'%kwargs['lam_drop'], '%.3f'%self.titer_tree.abs_error, '%.3f'%self.titer_tree.rms_error, '%.3f'%self.titer_tree.r2])+'\n')
            logfile.close()


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

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-v', '--viruses_per_month', type = int, default = 3, help='number of viruses sampled per month')
    parser.add_argument('-t', '--time_interval', nargs=2, default = ('1900-01-01', '2017-01-01'), help='specify time interval rather than use --years_back')
    parser.add_argument('-s', '--serotype', type = str, choices=['1', '2', '3', '4', 'any', 'all'], default='any', help = 'which serotype of dengue to build trees for; any = include all serotypes in one build; all = run all 5 builds')
    parser.add_argument('--load', action='store_true', help = 'recover from file')
    parser.add_argument('--fname', type = str, default=None, help = 'Specify input data filename (instead of inferring from serotype). Should be located in ../fauna/data directory.')
    parser.add_argument('--path_name', type = str, default=None, help = 'Prefix for output directories')
    parser.add_argument('--lam_pot', type = float, default=0.5, help = 'Regularization term: serum potency, Pb (L2)')
    parser.add_argument('--lam_avi', type = float, default=3.0, help = 'Regularization term: virus avidity, Va (L2)')
    parser.add_argument('--lam_drop', type = float, default=1.0, help = 'Regularization term: distance, Dab (L1)')
    parser.add_argument('--training_fraction', type = float, default=1.0, help = 'Fraction of titer measurements to use as training set.')

    params = parser.parse_args()

    # For any-serotype build, use dengue 4 outgroup and look for files like dengue.fasta
    serotype = params.serotype

    if serotype == 'any':
        lineage = 'dengue'
        reference_fname = './dengue/metadata/dengue_%s_outgroup.gb'%'4'
        sequence_fname = sorted(glob('../fauna/data/%s.fasta'%lineage), key=lambda f: os.path.getmtime(f))[-1]

    # For serotype-specific build, use the corresponding outgroup
    else:
        lineage = 'dengue_%s'%serotype
        reference_fname = './dengue/metadata/dengue_%s_outgroup.gb'%serotype
        try:    # Look for a serotype-specific fasta
            sequence_fname = sorted(glob('../fauna/data/%s*.fasta'%lineage), key=lambda f: os.path.getmtime(f))[-1]
        except: # If it doesn't exist, try to pull serotype-specific sequences out of the any-serotype fasta (warn the user of this behavior)
            sequence_fname = select_serotype('../fauna/data/dengue.fasta', '../fauna/data/', serotype)
            print('WARNING: Did not find serotype-specific fasta file.\nPulled sequences with serotype %s from any-serotype fasta file %s\nWrote these to file %s'%(serotype, '../fauna/data/dengue.fasta', newest_sequence_file))

    if params.fname:
        sequence_fname = '../fauna/data/'+params.fname

    input_data_path = sequence_fname.split('.fasta')[0]
    titer_fname = '../fauna/data/dengue_titers.tsv'
    strain_fname = '../fauna/data/dengue_strains.tsv'

    if params.path_name:
        store_data_path = 'store/'+params.path_name+'_'
        build_data_path = 'build/'+params.path_name+'_'
    else:
        store_data_path = 'store/'+lineage + '_'
        build_data_path = 'build/'+lineage + '_'

    dengue = dengue_process(input_data_path = input_data_path, store_data_path = store_data_path,
                   build_data_path = build_data_path, reference=reference_fname,
                   proteins= ['C', 'M', 'E', 'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', '2K', 'NS4B', 'NS5'],
                   method='SLSQP', inertia=np.exp(-1.0/12.0), stiffness=2.*12.0, titer_fname = titer_fname, strain_fname = strain_fname)

    if params.load: # Load sequences and tree from file
        dengue.load()
    else:
        dengue.load_sequences(fields={0:'strain', 1:'accession', 2:'date', 3:'region', 4:'country',
                        5:'division', 6: 'location'})

        time_interval = [datetime.strptime(x, '%Y-%m-%d').date() for x in params.time_interval]
        dropped_strains = ['DENV1/VIETNAM/BIDV992/2006', 'DENV1/FRANCE/00475/2008', 'DENV1/VIETNAM/BIDV3990/2008', 'DENV2/HAITI/DENGUEVIRUS2HOMOSAPIENS1/2016',# probable recombinants
        'DENV2/AUSTRALIA/QML22/2015' # Wrong serotype? Abnormal date.
        ]

        dengue.seqs.filter(lambda s: s.id not in dropped_strains)
        dengue.seqs.filter(lambda s:
            s.attributes['date']>=time_interval[0] and s.attributes['date']<time_interval[1])
        dengue.seqs.filter(lambda s: len(s.seq)>=900)
        dengue.subsample(params.viruses_per_month)
        dengue.align(debug=False)
        dengue.dump()

        dengue.align(debug=True)
        dengue.build_tree(debug=True)
        dengue.clock_filter(n_iqd=3, plot=True)
        dengue.annotate_tree(Tc=False, timetree=True, reroot='best')
        dengue.dump()

    dengue.titer_model(lam_pot = params.lam_pot, lam_avi = params.lam_avi, lam_drop = params.lam_drop, training_fraction = params.training_fraction)
    dengue.dump()
    dengue.titer_export()
    dengue.export(extra_attr=['serum'], controls=attribute_nesting, color_options=color_options, panels=panels )
