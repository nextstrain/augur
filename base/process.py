from __future__ import division, print_function
import sys, os, time, gzip, glob
from collections import defaultdict
from base.io_util import make_dir, remove_dir, tree_to_json, write_json, myopen
from base.sequences import sequence_set, num_date
from base.tree import tree
from base.frequencies import alignment_frequencies, tree_frequencies, make_pivots
import numpy as np
from datetime import datetime

class process(object):
    """process influenza virus sequences in mutliple steps to allow visualization in browser
        * filtering and parsing of sequences
        * alignment
        * tree building
        * frequency estimation of clades and mutations
        * export as json
    """

    def __init__(self, input_data_path = 'data/test',
                 store_data_path = 'store/test',
                 build_data_path = 'build/test', **kwargs):
        super(process, self).__init__()
        print("Initializing process")
        for p in [input_data_path, store_data_path, build_data_path]:
            if not os.path.isdir(os.path.dirname(p)):
                os.makedirs(os.path.dirname(p))
        self.input_data_path = input_data_path
        self.store_data_path = store_data_path
        self.build_data_path = build_data_path
        self.kwargs = kwargs
        self.data_filenames()

        self.seqs=None

        # parse the outgroup information
        if 'reference' in kwargs:
            reference_file = kwargs['reference']
            self.load_reference(reference_file)
        else:
            self.reference_seq = None


    def load_reference(self, reference_file):
        from Bio import SeqIO
        from Bio.SeqFeature import FeatureLocation
        self.reference_seq = SeqIO.read(reference_file, 'genbank')
        self.reference_seq.id = self.reference_seq.name
        self.genome_annotation = self.reference_seq.features
        if "proteins" in self.kwargs:
            # grap annotation from genbank
            protein_list = self.kwargs['proteins']
            self.proteins = {f.qualifiers['gene'][0]:FeatureLocation(start=f.location.start, end=f.location.end, strand=1)
                            for f in self.genome_annotation
                                if 'gene' in f.qualifiers
                                    and f.qualifiers['gene'][0] in protein_list}
        else:
            self.proteins = {}


    def load_sequences(self, fields={0:'strain', 2:'isolate_id', 3:'date', 4:'region',
                                     5:'country', 7:"city", 12:"subtype",13:'lineage'}):
        # instantiate and population the sequence objects
        self.seqs = sequence_set(self.sequence_fname, reference_seq=self.reference_seq)
        self.seqs.ungap()
        self.seqs.parse(fields, strip='_')

        # make sure the reference is part of the sequence set
        if self.reference_seq is not None:
            if self.reference_seq.name in self.seqs.all_seqs:
                self.seqs.all_seqs[self.reference_seq.name].seq=self.reference_seq.seq
            else:
                print('Outgroup is not in data base')

            # throw out sequences without dates
        self.seqs.parse_date(["%Y-%m-%d"], prune=True)


    def data_filenames(self):
        '''
        define filenames of input files and intermediates outputs
        '''
        self.sequence_fname = self.input_data_path+'.fasta'
        self.file_dumps = {}
        self.file_dumps['seqs'] = self.store_data_path+'sequences.pkl.gz'
        self.file_dumps['tree'] = self.store_data_path+'tree.newick'
        self.file_dumps['nodes'] = self.store_data_path+'nodes.pkl.gz'


    def dump(self):
        '''
        write the current state to file
        '''
        from cPickle import dump
        from Bio import Phylo
        for attr_name, fname in self.file_dumps.iteritems():
            if hasattr(self,attr_name):
                print("dumping",attr_name)
                #if attr_name=='seqs': self.seqs.all_seqs = None
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


    def align(self, outgroup=None):
        '''
        align sequences, remove non-reference insertions, outlier sequences, and translate
        '''
        self.seqs.align()
        self.seqs.strip_non_reference()
        if outgroup is not None:
            self.seqs.clock_filter(n_iqd=3, plot=False, max_gaps=0.05, root_seq=outgroup)
        self.seqs.translate(proteins=self.proteins)


    def estimate_mutation_frequencies(self, region="global", pivots=24):
        '''
        calculate the frequencies of mutation in a particular region
        currently the global frequencies should be estimated first
        because this defines the set of positions at which frequencies in
        other regions are estimated.
        '''
        if not hasattr(self.seqs, 'aln'):
            print("Align sequences first")
            return
        def filter_alignment(aln, region=None, lower_tp=None, upper_tp=None):
            from Bio.Align import MultipleSeqAlignment
            tmp = aln
            if region is not None:
                if type(region)==str:
                    tmp = [s for s in tmp if s.attributes['region']==region]
                elif type(region)==list:
                    tmp = [s for s in tmp if s.attributes['region'] in region]
                else:
                    print("region must be string or list")
                    return
            if lower_tp is not None:
                tmp = [s for s in tmp if np.mean(s.attributes['num_date'])>=lower_tp]
            if upper_tp is not None:
                tmp = [s for s in tmp if np.mean(s.attributes['num_date'])<upper_tp]
            return MultipleSeqAlignment(tmp)

        if not hasattr(self, 'pivots'):
            tps = np.array([np.mean(x.attributes['num_date']) for x in self.seqs.seqs.values()])
            self.pivots=make_pivots(pivots, tps)
        else:
            print('estimate_mutation_frequencies: using self.pivots')

        if not hasattr(self, 'mutation_frequencies'):
            self.mutation_frequencies = {}
            self.mutation_frequency_confidence = {}
            self.mutation_frequency_counts = {}

        # loop over nucleotide sequences and translations and calcuate
        # region specific frequencies of mutations above a certain threshold
        if type(region)==str:
            region_name = region
            region_match = region
        elif type(region)==tuple:
            region_name=region[0]
            region_match=region[1]
        else:
            print ("region must be string or tuple")
            return
        for prot, aln in [('nuc',self.seqs.aln)]+ self.seqs.translations.items():
            if region_match=="global":
                tmp_aln = filter_alignment(aln, lower_tp=self.pivots[0], upper_tp=self.pivots[-1])
                include_set=[]
            else:
                tmp_aln = filter_alignment(aln, region=region_match, lower_tp=self.pivots[0], upper_tp=self.pivots[-1])
                include_set = set([pos for (pos, mut) in self.mutation_frequencies[('global', prot)]])
            time_points = [np.mean(x.attributes['num_date']) for x in tmp_aln]
            if len(time_points)==0:
                print('no samples in region', region_name, prot)
                self.mutation_frequency_counts[region_name]=np.zeros_like(self.pivots)
                continue

            aln_frequencies = alignment_frequencies(tmp_aln, time_points, self.pivots,
                                            ws=max(2,len(time_points)//10),
                                            **self.kwargs)
            aln_frequencies.mutation_frequencies(min_freq=0.01)
            self.mutation_frequencies[(region_name,prot)] = aln_frequencies.frequencies
            self.mutation_frequency_confidence[(region_name,prot)] = aln_frequencies.calc_confidence()
            self.mutation_frequency_counts[region_name]=aln_frequencies.counts


    def estimate_tree_frequencies(self, region='global', pivots=24):
        '''
        estimate frequencies of clades in the tree, possibly region specific
        '''
        if region=='global':
            node_filter_func = None
        else:
            node_filter_func = lambda x:x.attr['region']==region

        if not hasattr(self, 'pivots'):
            tps = np.array([x.attributes['num_date'] for x in self.seqs.seqs.values()])
            self.pivots=make_pivots(pivots, tps)
        else:
            print('estimate_tree_frequencies: using self.pivots')
        if not hasattr(self, 'tree_frequencies'):
            self.tree_frequencies = {}
            self.tree_frequency_confidence = {}
            self.tree_frequency_counts = {}

        tree_freqs = tree_frequencies(self.tree.tree, self.pivots,
                                      node_filter = node_filter_func,
                                      ws = max(2,self.tree.tree.count_terminals()//10),
                                      **self.kwargs)

        tree_freqs.estimate_clade_frequencies()
        conf = tree_freqs.calc_confidence()
        self.tree_frequencies[region] = tree_freqs.frequencies
        self.tree_frequency_confidence[region] = conf
        self.tree_frequency_counts[region] = tree_freqs.counts


    def build_tree(self, infile=None, nodefile=None, root='best'):
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


    def clock_filter(self, n_iqd=3, plot=True):
        self.tree.tt.clock_filter(reroot='best', n_iqd=n_iqd, plot=plot)
        leaves = [x for x in self.tree.tree.get_terminals()]
        for n in leaves:
            if n.bad_branch:
                self.tree.tt.tree.prune(n)
                print('pruning leaf ', n.name)
        self.tree.tt.prepare_tree()


    def matchClades(self, clades, offset=-1):
        '''
        finds branches in the tree corresponding to named clades by searching for the
        oldest node with a particular genotype.
        - params
            - clades: a dictionary with clade names as keys and lists of genoypes as values
            - offset: the offset to be applied to the position specification, typically -1
                      to conform with counting starting at 0 as opposed to 1
        '''
        def match(node, genotype):
            return all([node.translations[gene][pos+offset]==state if gene in node.translations else node.sequences[pos+offset]==state
                        for gene, pos, state in genotype])

        self.clades_to_nodes = {}
        for clade_name, genotype in clades.iteritems():
            matching_nodes = filter(lambda x:match(x,genotype), self.tree.tree.get_nonterminals())
            matching_nodes.sort(key=lambda x:x.numdate if hasattr(x,'numdate') else x.dist2root)
            if len(matching_nodes):
                self.clades_to_nodes[clade_name] = matching_nodes[0]
                self.clades_to_nodes[clade_name].attr['clade_name']=clade_name
            else:
                print('matchClades: no match found for ', clade_name, genotype)


    def annotate_tree(self, Tc=0.01, timetree=False, **kwargs):
        if timetree:
            self.tree.timetree(Tc=Tc, infer_gtr=True, **kwargs)
        else:
            self.tree.ancestral(**kwargs)
        self.tree.add_translations()
        self.tree.refine()
        self.tree.layout()


    def make_control_json(self, controls):
        controls_json = {}
        for super_cat, fields in controls.iteritems():
            cat_count = {}
            for n in self.tree.tree.get_terminals():
                tmp = cat_count
                for field in fields:
                    tmp["name"] = field
                    if field in n.attr:
                        cat = n.attr[field]
                    else:
                        cat='unknown'
                    if cat in tmp:
                        tmp[cat]['count']+=1
                    else:
                        tmp[cat] = {'count':1, 'subcats':{}}
                    tmp = tmp[cat]['subcats']
            controls_json[super_cat] = cat_count
        return controls_json


    def export(self, extra_attr = [], controls = {}):
        '''
        export the tree, sequences, frequencies to json files for visualization
        in the browser
        '''
        prefix = self.build_data_path
        # export json file that contains alignment diversity column by column
        self.seqs.export_diversity(prefix+'entropy.json')
        # exports the tree and the sequences inferred for all clades in the tree
        if hasattr(self, 'tree') and self.tree is not None:
            self.tree.export(path=prefix, extra_attr = extra_attr
                         + ["muts", "aa_muts","attr", "clade"])


        # local function or round frequency estimates to useful precision (reduces file size)
        def process_freqs(freq):
            return [round(x,4) for x in freq]

        # construct a json file containing all frequency estimate
        # the format is region_protein:159F for mutations and region_clade:123 for clades
        if hasattr(self, 'pivots'):
            freq_json = {'pivots':process_freqs(self.pivots)}
        if hasattr(self, 'mutation_frequencies'):
            freq_json['counts'] = {x:list(counts) for x, counts in self.mutation_frequency_counts.iteritems()}
            for (region, gene), tmp_freqs in self.mutation_frequencies.iteritems():
                for mut, freq in tmp_freqs.iteritems():
                    label_str =  region+"_"+ gene + ':' + str(mut[0]+1)+mut[1]
                    freq_json[label_str] = process_freqs(freq)
        # repeat for clade frequencies in trees
        if hasattr(self, 'tree_frequencies'):
            for region in self.tree_frequencies:
                for clade, freq in self.tree_frequencies[region].iteritems():
                    label_str = region+'_clade:'+str(clade)
                    freq_json[label_str] = process_freqs(freq)
        # repeat for named clades
        if hasattr(self, 'clades_to_nodes') and hasattr(self, 'tree_frequencies'):
            for region in self.tree_frequencies:
                for clade, node in self.clades_to_nodes.iteritems():
                    label_str = region+'_'+str(clade)
                    freq_json[label_str] = process_freqs(self.tree_frequencies[region][node.clade])
        # write to one frequency json
        if hasattr(self, 'tree_frequencies') or hasattr(self, 'mutation_frequencies'):
            write_json(freq_json, prefix+'frequencies.json', indent=None)
        if len(controls):
            controls_json = self.make_control_json(controls)
            write_json(controls_json, prefix+'controls.json')

        # write out metadata json# Write out metadata
        print("Writing out metadata")
        meta_json = {}
        meta_json["updated"] = time.strftime("X%d %b %Y").replace('X0','X').replace('X','')
        try:
            from pygit2 import Repository, discover_repository
            current_working_directory = os.getcwd()
            repository_path = discover_repository(current_working_directory)
            repo = Repository(repository_path)
            commit_id = repo[repo.head.target].id
            meta_json["commit"] = str(commit_id)
        except ImportError:
            meta_json["commit"] = "unknown"
        write_json(meta_json, prefix+'meta.json')


if __name__=="__main__":
    lineage = 'h3n2'
    input_data_path = '../nextstrain-db/data/'+lineage
    store_data_path = 'store/'+lineage + '_'
    build_data_path = 'build/'+lineage + '_'

    proc = process(input_data_path = input_data_path, store_data_path = store_data_path, build_data_path = build_data_path,
                   reference='flu/metadata/h3n2_outgroup.gb', proteins=['HA1', 'HA2'],method='SLSQP')
    proc.load_sequences()
    proc.seqs.filter(lambda s: s.attributes['date']>=datetime(2012,1,1).date() and
                               s.attributes['date']< datetime(2016,1,1).date())
    proc.seqs.subsample(category = lambda x:(x.attributes['region'],
                                             x.attributes['date'].year,
                                             x.attributes['date'].month), threshold=1)

    proc.align()
    proc.estimate_mutation_frequencies(region='global')
    proc.build_tree()
    proc.annotate_tree(Tc=0.005, timetree=True)
    proc.estimate_tree_frequencies(region='global')
    proc.estimate_tree_frequencies(region='north_america')
    proc.tree.geo_inference('region')
    proc.export()
