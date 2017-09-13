from __future__ import division, print_function
import argparse
import sys, os, time, gzip, glob
from collections import defaultdict
from base.config import combine_configs
from base.io_util import make_dir, remove_dir, tree_to_json, write_json, myopen
from base.sequences_process import sequence_set
from base.utils import num_date, save_as_nexus, parse_date
from base.tree import tree
from base.fitness_model import fitness_model
from base.frequencies import alignment_frequencies, tree_frequencies, make_pivots
from base.auspice_export import export_metadata_json, export_frequency_json
import numpy as np
from datetime import datetime
import json
from pdb import set_trace
from base.logger import logger
from Bio import SeqIO
from Bio import AlignIO
import cPickle as pickle


def collect_args():
    parser = argparse.ArgumentParser(
        description = "Process (prepared) JSON(s)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-j', '--json', help="prepared JSON to process")
    parser.add_argument('--clean', default=False, action='store_true', help="clean build (remove previous checkpoints)")
    parser.add_argument('--no_raxml', action='store_true', help="do not run RAxML to build the tree")
    parser.add_argument('--no_tree', action='store_true', help="do not build a tree")

    return parser


class process(object):
    """process influenza virus sequences in mutliple steps to allow visualization in browser
        * filtering and parsing of sequences
        * alignment
        * tree building
        * frequency estimation of clades and mutations
        * export as json
    """

    def __init__(self, config):
        """ check config file, make necessary directories, set up logger """
        super(process, self).__init__()
        self.config = combine_configs("process", config)

        try:
            assert(os.path.basename(os.getcwd()) == self.config["dir"])
        except AssertionError:
            print("Run this script from within the {} directory".format(self.config["dir"]))
            sys.exit(2)

        for p in self.config["output"].values():
            if not os.path.isdir(p):
                os.makedirs(p)

        self.log = logger(self.config["output"]["data"], False)

        # parse the JSON into different data bits
        try:
            with open(self.config["in"], 'r') as fh:
                data = json.load(fh)
        except Exception as e:
            self.log.fatal("Error loading JSON. Error: {}".format(e))

        self.info = data["info"]
        if "time_interval" in data["info"]:
            self.info["time_interval"] = [datetime.strptime(x, '%Y-%m-%d').date()
                                          for x in data["info"]["time_interval"]]
        self.info["lineage"] = data["info"]["lineage"]

        if 'leaves' in data:
            self.tree_leaves = data['leaves']

        try:
            self.colors = data["colors"]
        except KeyError:
            self.log.notify("* colours have not been set")
            self.colors = False
        try:
            self.lat_longs = data["lat_longs"]
        except KeyError:
            self.log.notify("* latitude & longitudes have not been set")
            self.lat_longs = False

        # backwards compatability - set up file_dumps (need to rewrite sometime)
        # self.sequence_fname = self.input_data_path+'.fasta'
        self.file_dumps = {}
        self.output_path = os.path.join(self.config["output"]["data"], self.info["prefix"])
        self.file_dumps['seqs'] = self.output_path + '_sequences.pkl.gz'
        self.file_dumps['tree'] = self.output_path + '_tree.newick'
        self.file_dumps['nodes'] = self.output_path + '_nodes.pkl.gz'

        if self.config["clean"] == True:
            self.log.notify("Removing intermediate files for a clean build")
            for f in glob.glob(self.output_path+"*"):
                os.remove(f)

        if "reference" in data:
            self.seqs = sequence_set(self.log, data["sequences"], data["reference"], self.info["date_format"])
        else:
            self.log.fatal("No reference provided. Cannot continue.")
            # self.seqs = sequence_set(self.log, data["sequences"], False, self.info["date_format"])
        # backward compatability
        self.reference_seq = self.seqs.reference_seq
        self.proteins = self.seqs.proteins

        for trait in self.info["traits_are_dates"]:
            self.seqs.convert_trait_to_numerical_date(trait, self.info["date_format"])

        # Prepare titers if they are available.
        if "titers" in data:
            self.log.debug("Loaded %i titer measurements" % len(data["titers"]))
            # Convert titer dictionary indices from JSON-compatible strings back
            # to tuples.
            self.titers = {eval(key): value
                           for key, value in data["titers"].iteritems()}

        ## usefull flag to set (from pathogen run file) to disable restoring
        self.try_to_restore = True


    def dump(self):
        '''
        write the current state to file
        '''
        self.log.warn("unsure if dump() works")
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

    def load(self, debug=False):
        '''
        reconstruct instance from files
        '''
        self.log.warn("unsure if load() works")
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
            self.build_tree(tree_name, node_file, root='none', debug=debug)


    def align(self, codon_align=False, debug=False, fill_gaps=False):
        '''
        (1) Align sequences, remove non-reference insertions
        NB step 1 is skipped if a valid aln file is found
        (2) Translate
        (3) Write to multi-fasta
        CODON ALIGNMENT IS NOT IMPLEMENTED
        '''
        fnameStripped = self.output_path + "_aligned_stripped.mfa"
        if self.try_to_restore:
            self.seqs.try_restore_align_from_disk(fnameStripped)
        if not hasattr(self.seqs, "aln"):
            if codon_align:
                self.seqs.codon_align()
            else:
                self.seqs.align(self.config["subprocess_verbosity_level"], debug=debug)
            # need to redo everything
            self.try_to_restore = False

            self.seqs.strip_non_reference()
            if fill_gaps:
                self.seqs.make_gaps_ambiguous()
            else:
                self.seqs.make_terminal_gaps_ambiguous()


            AlignIO.write(self.seqs.aln, fnameStripped, 'fasta')

            if not self.seqs.reference_in_dataset:
                self.seqs.remove_reference_from_alignment()
            # if outgroup is not None:
            #     self.seqs.clock_filter(n_iqd=3, plot=False, max_gaps=0.05, root_seq=outgroup)

        self.seqs.translate() # creates self.seqs.translations
        # save additional translations - disabled for now
        # for name, msa in self.seqs.translations.iteritems():
        #     SeqIO.write(msa, self.output_path + "_aligned_" + name + ".mfa", "fasta")


    def get_pivots_via_spacing(self):
        try:
            time_interval = self.info["time_interval"]
            assert("pivot_spacing" in self.config)
        except AssertionError:
            self.log.fatal("Cannot space pivots without prividing \"pivot_spacing\" in the config")
        except KeyError:
            self.log.fatal("Cannot space pivots without a time interval in the prepared JSON")
        return np.arange(time_interval[1].year+(time_interval[1].month-1)/12.0,
                         time_interval[0].year+time_interval[0].month/12.0,
                         self.config["pivot_spacing"])

    def restore_mutation_frequencies(self):
        if self.try_to_restore:
            try:
                with open(self.output_path + "_mut_freqs.pickle", 'rb') as fh:
                    pickle_seqs = pickle.load(fh)
                    assert(pickle_seqs == set(self.seqs.seqs.keys()))
                    pickled = pickle.load(fh)
                    assert(len(pickled) == 3)
                    self.mutation_frequencies = pickled[0]
                    self.mutation_frequency_confidence = pickled[1]
                    self.mutation_frequency_counts = pickled[2]
                    self.log.notify("Successfully restored mutation frequencies")
                    return
            except IOError:
                pass
            except AssertionError as err:
                self.log.notify("Tried to restore mutation frequencies but failed: {}".format(err))
            #no need to remove - we'll overwrite it shortly
        self.mutation_frequencies = {}
        self.mutation_frequency_confidence = {}
        self.mutation_frequency_counts = {}

    def estimate_mutation_frequencies(self,
                                      inertia=0.0,
                                      min_freq=0.01,
                                      stiffness=20.0,
                                      pivots=24,
                                      region="global",
                                      include_set={}):
        '''
        calculate the frequencies of mutation in a particular region
        currently the global frequencies should be estimated first
        because this defines the set of positions at which frequencies in
        other regions are estimated.
        '''
        if not hasattr(self.seqs, 'aln'):
            self.log.warn("Align sequences first")
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
                    self.log.warn("region must be string or list")
                    return
            if lower_tp is not None:
                tmp = [s for s in tmp if np.mean(s.attributes['num_date'])>=lower_tp]
            if upper_tp is not None:
                tmp = [s for s in tmp if np.mean(s.attributes['num_date'])<upper_tp]
            return MultipleSeqAlignment(tmp)

        if not hasattr(self, 'pivots'):
            tps = np.array([np.mean(x.attributes['num_date']) for x in self.seqs.seqs.values()])
            self.pivots=make_pivots(pivots, tps)
        # else:
        #     self.log.notify('estimate_mutation_frequencies: using self.pivots')

        if not hasattr(self, 'mutation_frequencies'):
            self.restore_mutation_frequencies()

        # loop over nucleotide sequences and translations and calcuate
        # region specific frequencies of mutations above a certain threshold
        if type(region)==str:
            region_name = region
            region_match = region
        elif type(region)==tuple:
            region_name=region[0]
            region_match=region[1]
        else:
            self.log.warn("region must be string or tuple")
            return

        # loop over different alignment types
        for prot, aln in [('nuc',self.seqs.aln)] + self.seqs.translations.items():
            if (region_name,prot) in self.mutation_frequencies:
                self.log.notify("Skipping Frequency Estimation for region \"{}\", protein \"{}\"".format(region_name, prot))
                continue
            self.log.notify("Starting Frequency Estimation for region \"{}\", protein \"{}\"".format(region_name, prot))

            # determine set of positions that have to have a frequency calculated
            if prot in include_set:
                tmp_include_set = [x for x in include_set[prot]]
            else:
                tmp_include_set = []

            tmp_aln = filter_alignment(aln, region = None if region=='global' else region_match,
                                      lower_tp=self.pivots[0], upper_tp=self.pivots[-1])

            if ('global', prot) in self.mutation_frequencies:
                tmp_include_set += set([pos for (pos, mut) in self.mutation_frequencies[('global', prot)]])

            time_points = [np.mean(x.attributes['num_date']) for x in tmp_aln]
            if len(time_points)==0:
                self.log.notify('no samples in region {} (protein: {})'.format(region_name, prot))
                self.mutation_frequency_counts[region_name]=np.zeros_like(self.pivots)
                continue

            # instantiate alignment frequency
            aln_frequencies = alignment_frequencies(tmp_aln, time_points, self.pivots,
                                            ws=max(2,len(time_points)//10),
                                            inertia=inertia,
                                            stiffness=stiffness, method='SLSQP')
            if prot=='nuc': # if this is a nucleotide alignment, set all non-canonical states to N
                A = aln_frequencies.aln
                A[~((A=='A')|(A=='C')|(A=='G')|(A=='T')|('A'=='-'))] = 'N'

            aln_frequencies.mutation_frequencies(min_freq=min_freq, include_set=tmp_include_set,
                                                 ignore_char='N' if prot=='nuc' else 'X')
            self.mutation_frequencies[(region_name,prot)] = aln_frequencies.frequencies
            self.mutation_frequency_confidence[(region_name,prot)] = aln_frequencies.calc_confidence()
            self.mutation_frequency_counts[region_name]=aln_frequencies.counts

        self.log.notify("Saving mutation frequencies (pickle)")
        with open(self.output_path + "_mut_freqs.pickle", 'wb') as fh:
            pickle.dump(set(self.seqs.seqs.keys()), fh, protocol=pickle.HIGHEST_PROTOCOL)
            pickle.dump((self.mutation_frequencies,
                         self.mutation_frequency_confidence,
                         self.mutation_frequency_counts), fh, protocol=pickle.HIGHEST_PROTOCOL)


    def global_frequencies(self, min_freq):
        # determine sites whose frequencies need to be computed in all regions
        self.seqs.diversity_statistics()
        include_set = {}
        for prot in ['nuc'] + self.seqs.translations.keys():
            include_set[prot] = np.where(np.sum(self.seqs.af[prot][:-2]**2, axis=0)<np.sum(self.seqs.af[prot][:-2], axis=0)**2-min_freq)[0]

        # set pivots and define groups of larger regions for frequency display
        pivots = self.get_pivots_via_spacing()
        acronyms = set([x[1] for x in self.info["regions"] if x[1]!=""])
        region_groups = {str(x):[str(y[0]) for y in self.info["regions"] if y[1] == x] for x in acronyms}
        pop_sizes = {str(x):np.sum([y[-1] for y in self.info["regions"] if y[1] == x]) for x in acronyms}
        total_popsize = np.sum(pop_sizes.values())

        # estimate frequencies in individual regions
        # TODO: move inertia and stiffness parameters to config
        for region in region_groups.iteritems():
            self.estimate_mutation_frequencies(pivots=pivots, region=region, min_freq=0.02, include_set=include_set,
                                                 inertia=np.exp(-2.0/12), stiffness=2.0*12)

        # perform a weighted average of frequencies across the regions to determine
        # global frequencies.
        # First: compute the weights accounting for seasonal variation and populations size
        weights = {region: np.array(self.mutation_frequency_counts[region], dtype = float)
                   for region in acronyms}

        for region in weights: # map maximal count across time to 1.0, weigh by pop size
            weights[region] = np.maximum(0.1, weights[region]/weights[region].max())
            weights[region]*=pop_sizes[region]

        # compute the normalizer
        total_weight = np.sum([weights[region] for region in acronyms],axis=0)

        for prot in ['nuc'] + self.seqs.translations.keys():
            gl_freqs, gl_counts, gl_confidence = {}, {}, {}
            all_muts = set()
            for region in acronyms: # list all unique mutations
                all_muts.update(self.mutation_frequencies[(region, prot)].keys())
            for mut in all_muts: # compute the weighted average
                gl_freqs[mut] = np.sum([self.mutation_frequencies[(region, prot)][mut]*weights[region] for region in acronyms
                                        if mut in self.mutation_frequencies[(region, prot)]], axis=0)/total_weight
                gl_confidence[mut] = np.sqrt(np.sum([self.mutation_frequency_confidence[(region, prot)][mut]**2*weights[region]
                                                     for region in acronyms
                                            if mut in self.mutation_frequencies[(region, prot)]], axis=0)/total_weight)
            gl_counts = np.sum([self.mutation_frequency_counts[region] for region in acronyms
                                        if mut in self.mutation_frequencies[(region, prot)]], axis=0)
            # save in mutation_frequency data structure
            self.mutation_frequencies[("global", prot)] = gl_freqs
            self.mutation_frequency_counts["global"] = gl_counts
            self.mutation_frequency_confidence[("global", prot)] = gl_confidence


    def save_tree_frequencies(self):
        """
        Save tree frequencies to a pickle on disk.
        """
        self.log.notify("Saving tree frequencies (pickle)")
        with open(self.output_path + "_tree_freqs.pickle", 'wb') as fh:
            pickle.dump(set(self.seqs.seqs.keys()), fh, protocol=pickle.HIGHEST_PROTOCOL)
            pickle.dump((self.tree_frequencies,
                         self.tree_frequency_confidence,
                         self.tree_frequency_counts,
                         self.pivots), fh, protocol=pickle.HIGHEST_PROTOCOL)


    def restore_tree_frequencies(self):
        try:
            assert(self.try_to_restore == True)
            with open(self.output_path + "_tree_freqs.pickle", 'rb') as fh:
                pickle_seqs = pickle.load(fh)
                assert(pickle_seqs == set(self.seqs.seqs.keys()))
                pickled = pickle.load(fh)
                assert(len(pickled) == 4)
                self.tree_frequencies = pickled[0]
                self.tree_frequency_confidence = pickled[1]
                self.tree_frequency_counts = pickled[2]
                self.pivots = pickled[3]
                self.log.notify("Successfully restored tree frequencies")
                return
        except IOError:
            pass
        except AssertionError as err:
            self.log.notify("Tried to restore tree frequencies but failed: {}".format(err))
            #no need to remove - we'll overwrite it shortly
        self.tree_frequencies = {}
        self.tree_frequency_confidence = {}
        self.tree_frequency_counts = {}


    def estimate_tree_frequencies(self, region='global', pivots=24):
        '''
        estimate frequencies of clades in the tree, possibly region specific
        '''
        if region=='global':
            node_filter_func = None
        else:
            node_filter_func = lambda x:x.attr['region']==region

        if not hasattr(self, 'tree_frequencies'):
            self.restore_tree_frequencies()

        if region in self.tree_frequencies:
            self.log.notify("Skipping tree frequency estimation for region: %s" % region)
            return

        if not hasattr(self, 'pivots'):
            tps = np.array([np.mean(x.attributes['num_date']) for x in self.seqs.seqs.values()])
            self.pivots=make_pivots(pivots, tps)

        self.log.notify('Estimate tree frequencies for %s: using self.pivots' % (region))

        tree_freqs = tree_frequencies(self.tree.tree, self.pivots, method='SLSQP',
                                      node_filter = node_filter_func,
                                      ws = max(2,self.tree.tree.count_terminals()//10))
                                    # who knows what kwargs are needed here
                                    #   **self.kwargs)

        tree_freqs.estimate_clade_frequencies()
        conf = tree_freqs.calc_confidence()
        self.tree_frequencies[region] = tree_freqs.frequencies
        self.tree_frequency_confidence[region] = conf
        self.tree_frequency_counts[region] = tree_freqs.counts

        self.save_tree_frequencies()

    def build_tree(self):
        '''
        (1) instantiate a tree object (process.tree)
        (2) If newick file doesn't exist or isn't valid: build a newick tree (normally RAxML)
        (3) Make a TimeTree
        '''
        self.tree = tree(aln=self.seqs.aln, proteins=self.proteins, verbose=self.config["subprocess_verbosity_level"])
        newick_file = self.output_path + ".newick"
        if self.try_to_restore and os.path.isfile(newick_file) and self.tree.check_newick(newick_file):
            self.log.notify("Newick file \"{}\" can be used to restore".format(newick_file))
        else:
            self.log.notify("Building newick tree.")
            self.tree.build_newick(newick_file, **self.config["newick_tree_options"])


    def clock_filter(self):
        if self.config["clock_filter"] == False:
            return
        self.tree.tt.clock_filter(reroot='best', n_iqd=self.config["clock_filter"]["n_iqd"], plot=self.config["clock_filter"]["plot"])

        leaves = [x for x in self.tree.tree.get_terminals()]
        for n in leaves:
            if n.bad_branch:
                self.tree.tt.tree.prune(n)
                print('pruning leaf ', n.name)
        if self.config["clock_filter"]["remove_deep_splits"]:
            self.tree.tt.tree.ladderize()
            current_root = self.tree.tt.tree.root
            if sum([x.branch_length for x in current_root])>0.1 \
                and sum([x.count_terminals() for x in current_root.clades[:-1]])<5:
                new_root = current_root.clades[-1]
                new_root.up=False
                self.tree.tt.tree.root = new_root
                with open(self.output_path+"_outliers.txt", 'a') as ofile:
                    for x in current_root.clades[:-1]:
                        ofile.write("\n".join([leaf.name for leaf in x.get_terminals()])+'\n')

        self.tree.tt.prepare_tree()


    def timetree_setup_filter_run(self):
        def try_restore():
            try:
                assert(os.path.isfile(self.output_path + "_timetree.new"))
                assert(os.path.isfile(self.output_path + "_timetree.pickle"))
            except AssertionError:
                return False

            self.log.notify("Attempting to restore timetree")
            with open(self.output_path+"_timetree.pickle", 'rb') as fh:
                pickled = pickle.load(fh)

            try:
                assert(self.config["timetree_options"] == pickled["timetree_options"])
                assert(self.config["clock_filter"] == pickled["clock_filter_options"])
                #assert(set(self.seqs.sequence_lookup.keys()) == set(pickled["original_seqs"]))
            except AssertionError as e:
                print(e)
                self.log.warn("treetime is out of date - rerunning")
                return False

            # this (treetime) newick is _after_ clock filtering and remove_outliers_clades
            # so these methods should not be rerun here
            self.tree.tt_from_file(self.output_path + "_timetree.new", nodefile=None, root=None)
            try:
                self.tree.restore_timetree_node_info(pickled["nodes"])
            except KeyError:
                self.log.warn("treetime node info missing - rerunning")
                return False
            self.log.notify("TreeTime successfully restored.")
            return True

        if "temporal_confidence" in self.config:
            self.config["timetree_options"]["confidence"] = True
            self.config["timetree_options"]["use_marginal"] = True

        if self.try_to_restore:
            success = try_restore()
        else:
            success = False
        if not success:
            self.log.notify("Setting up TimeTree")
            self.tree.tt_from_file(self.output_path + ".newick", nodefile=None, root="best")
            self.log.notify("Running Clock Filter")
            self.clock_filter()
            self.tree.remove_outlier_clades() # this is deterministic
            self.log.notify("Reconstructing Ancestral Sequences, branch lengths & dating nodes")
            self.tree.timetree(**self.config["timetree_options"])
            # do we ever not want to use timetree?? If so:
            # self.tree.ancestral(**kwargs) instead of self.tree.timetree
            self.tree.save_timetree(fprefix=self.output_path, ttopts=self.config["timetree_options"], cfopts=self.config["clock_filter"])

        self.tree.add_translations()
        self.tree.refine()
        self.tree.layout()


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
            return all([node.translations[gene][pos+offset]==state if gene in node.translations else node.sequence[pos+offset]==state
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
                for allele in genotype:
                    partial_matches = filter(lambda x:match(x,[allele]), self.tree.tree.get_nonterminals())
                    print('Found %d partial matches for allele '%len(partial_matches), allele)

    def annotate_fitness(self):
        """Run the fitness prediction model and annotate the tree's nodes with fitness
        values. Returns the resulting fitness model instance.
        """
        kwargs = {
            "tree": self.tree.tree,
            "frequencies": self.tree_frequencies,
            "time_interval": self.info["time_interval"],
            "pivot_spacing": self.config["pivot_spacing"]
        }

        if "predictors" in self.config:
            kwargs["predictor_input"] = self.config["predictors"]

        if "epitope_mask" in self.config:
            kwargs["epitope_masks_fname"] = self.config["epitope_mask"]

        if "epitope_mask_version" in self.config:
            kwargs["epitope_mask_version"] = self.config["epitope_mask_version"]

        if self.config["subprocess_verbosity_level"] > 0:
            kwargs["verbose"] = 1

        model = fitness_model(**kwargs)
        model.predict()

        return model

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

    def auspice_export(self):
        '''
        export the tree, sequences, frequencies to json files for auspice visualization
        '''
        prefix = os.path.join(self.config["output"]["auspice"], self.info["prefix"])
        indent = 2

        ## ENTROPY (alignment diversity) ##
        self.seqs.export_diversity(fname=prefix+'_entropy.json', indent=indent)

        ## TREE (includes inferred states, mutations etc) ##
        if hasattr(self, 'tree') and self.tree is not None:
            self.tree.export(path=prefix, extra_attr = self.config["auspice"]["extra_attr"]
                         + ["muts", "aa_muts","attr", "clade"], indent = indent)

        ## FREQUENCIES ##
        export_frequency_json(self, prefix=prefix, indent=indent)

        ## METADATA ##
        export_metadata_json(self, prefix=prefix, indent=indent)


    def run_geo_inference(self):
        if self.config["geo_inference"] == False:
            self.log.notify("Not running geo inference")
            return
        try:
            kwargs = {"report_confidence": self.config["geo_inference_options"]["confidence"]}
        except KeyError:
            kwargs = {}

        ## try load pickle...
        try:
            assert(self.try_to_restore == True)
            with open(self.output_path + "_mugration.pickle", 'rb') as fh:
                options = pickle.load(fh)
                restored_data = pickle.load(fh)
            assert(options == self.config["geo_inference_options"])
            assert(set(restored_data.keys()) == set([x.name for x in self.tree.tree.find_clades()]))
        except IOError:
            restored_data = False
        except AssertionError as err:
            restored_data = False
            self.log.notify("Tried to restore mutation frequencies but failed: {}".format(err))

        # only run geo inference if lat + longs are defined.
        if not self.lat_longs or len(self.lat_longs)==0:
            self.log.notify("no geo inference - no specified lat/longs")
            return
        for geo_attr in self.config["geo_inference"]:
            try:
                self.tree.restore_geo_inference(restored_data, geo_attr, self.config["geo_inference_options"]["confidence"])
                self.log.notify("Restored geo inference for {}".format(geo_attr))
            except KeyError:
                try:
                    kwargs["root_state"] = self.config["geo_inference_options"]["root_state"][geo_attr]
                except KeyError:
                    pass
                self.log.notify("running geo inference for {} with parameters {}".format(geo_attr, kwargs))
                self.tree.geo_inference(geo_attr, **kwargs)

        # SAVE MUGRATION RESULTS:
        attrs = set(self.tree.mugration_attrs)
        try:
            data = {}
            for node in self.tree.tree.find_clades():
                assert(len(attrs - set(node.attr.keys()))==0)
                data[node.name] = {x:node.attr[x] for x in attrs}
        except AssertionError:
            self.log.warn("Error saving mugration data - will not be able to restore")
            return
        with open(self.output_path + "_mugration.pickle", 'wb') as fh:
            pickle.dump(self.config["geo_inference_options"], fh, protocol=pickle.HIGHEST_PROTOCOL)
            pickle.dump(data, fh, protocol=pickle.HIGHEST_PROTOCOL)
        self.log.notify("Saved mugration data (pickle)")

    def save_as_nexus(self):
        save_as_nexus(self.tree.tree, self.output_path + "_timeTree.nex")

if __name__=="__main__":
    print("This shouldn't be called as a script.")
