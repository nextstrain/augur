from __future__ import division, print_function
import sys, os, time, gzip, glob
from collections import defaultdict
from base.io_util import make_dir, remove_dir, tree_to_json, write_json, myopen
from base.sequences_process import sequence_set
from base.utils import num_date, save_as_nexus, parse_date
from base.tree import tree
from base.frequencies import alignment_frequencies, tree_frequencies, make_pivots
import numpy as np
from datetime import datetime
import json
from pdb import set_trace
from base.logger import logger
from Bio import SeqIO

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
        self.config = config

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

        if "reference" in data:
            self.seqs = sequence_set(self.log, data["sequences"], data["reference"], self.info["date_format"])
        else:
            self.seqs = sequence_set(self.log, data["sequences"], False, self.info["date_format"])
        # backward compatability
        self.reference_seq = self.seqs.reference_seq
        self.proteins = self.seqs.proteins

        for trait in self.info["traits_are_dates"]:
            self.seqs.convert_trait_to_numerical_date(trait, self.info["date_format"])

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


    def align(self, codon_align=False, debug=False):
        '''
        align sequences, remove non-reference insertions, outlier sequences, and translate
        '''
        if codon_align:
            self.seqs.codon_align(debug=debug)
        else:
            self.seqs.align(debug=debug)
        self.seqs.strip_non_reference()
        self.seqs.remove_terminal_gaps()
        # if outgroup is not None:
        #     self.seqs.clock_filter(n_iqd=3, plot=False, max_gaps=0.05, root_seq=outgroup)
        self.seqs.translate() # creates self.seqs.translations
        # save the final alignment!
        SeqIO.write(self.seqs.aln, self.output_path + "_aligned.mfa", "fasta")
        for name, msa in self.seqs.translations.iteritems():
            SeqIO.write(msa, self.output_path + "_aligned_" + name + ".mfa", "fasta")


    def estimate_mutation_frequencies(self, region="global", pivots=24, include_set=None, min_freq=0.01):
        '''
        calculate the frequencies of mutation in a particular region
        currently the global frequencies should be estimated first
        because this defines the set of positions at which frequencies in
        other regions are estimated.
        '''
        if include_set is None:
            include_set = {}
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
            if prot in include_set:
                tmp_include_set = [x for x in include_set[prot]]
            else:
                tmp_include_set = []
            if region_match=="global":
                tmp_aln = filter_alignment(aln, lower_tp=self.pivots[0], upper_tp=self.pivots[-1])
            else:
                tmp_aln = filter_alignment(aln, region=region_match, lower_tp=self.pivots[0], upper_tp=self.pivots[-1])
                tmp_include_set += set([pos for (pos, mut) in self.mutation_frequencies[('global', prot)]])
            time_points = [np.mean(x.attributes['num_date']) for x in tmp_aln]
            if len(time_points)==0:
                print('no samples in region', region_name, prot)
                self.mutation_frequency_counts[region_name]=np.zeros_like(self.pivots)
                continue

            aln_frequencies = alignment_frequencies(tmp_aln, time_points, self.pivots,
                                            ws=max(2,len(time_points)//10),
                                            **self.kwargs)
            aln_frequencies.mutation_frequencies(min_freq=min_freq, include_set=tmp_include_set,
                                                 ignore_char='N' if prot=='nuc' else 'X')
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
            print('estimate_tree_frequencies: using self.pivots', self.pivots)
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


    def build_tree(self, infile=None, nodefile=None, root='best', debug=False, num_distinct_starting_trees=1):
        '''
        instantiate a tree object and make a time tree
        if infiles are None, the tree is build from scratch. Otherwise
        the tree is loaded from file
        '''
        self.log.warn("self.verbose not set")
        self.tree = tree(aln=self.seqs.aln, proteins=self.proteins, verbose=2)
        if infile is None:
            self.tree.build(root=root, debug=debug, num_distinct_starting_trees=num_distinct_starting_trees)
        else:
            self.tree.tt_from_file(infile, nodefile=nodefile, root=root)


    def clock_filter(self, n_iqd=3, plot=True, remove_deep_splits=False):
        self.tree.tt.clock_filter(reroot='best', n_iqd=n_iqd, plot=plot)

        leaves = [x for x in self.tree.tree.get_terminals()]
        for n in leaves:
            if n.bad_branch:
                self.tree.tt.tree.prune(n)
                print('pruning leaf ', n.name)
        if remove_deep_splits:
            self.tree.tt.tree.ladderize()
            current_root = self.tree.tt.tree.root
            if sum([x.branch_length for x in current_root])>0.1 \
                and sum([x.count_terminals() for x in current_root.clades[:-1]])<5:
                new_root = current_root.clades[-1]
                new_root.up=False
                self.tree.tt.tree.root = new_root
                with open(self.store_data_path+"outliers.txt", 'a') as ofile:
                    for x in current_root.clades[:-1]:
                        ofile.write("\n".join([leaf.name for leaf in x.get_terminals()])+'\n')

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

    def define_latitude_longitude(self):
        import csv
        # get the latitude and longitudes that were already determined
        file = open(self.lat_long_fname, 'r')
        reader = csv.DictReader(filter(lambda row: row[0]!='#', file), delimiter='\t')		# list of dicts
        self.location_to_lat_long = {}
        for line in reader:
            try:
                self.location_to_lat_long[line['location']] = {
                    'latitude': float(line['latitude']),
                    'longitude': float(line['longitude'])
                    }
            except:
                print("Line failed ", line)
                raise Exception("Failed to read ", file, "please check the line that failed")
        file.close()

    # def make_geo_lookup_json(self, geo_attributes = []):
    #     '''
    #     Take existing geo attributes (region, country, division) for viruses and
    #     produces a lookup JSON to go from geo string to lat/long.
    #     Example:
    #     "geo_lookup": {
    #         "country": {
    #             "brazil": {
    #                 "latitude": -10.3333332,
    #                 "longitude": -53.1999999,
    #             },
    #             "colombia": {
    #                 "latitude": 2.893108,
    #                 "longitude": -73.7845142,
    #             }
    #         },
    #         "region": {
    #             "north_america": {
    #                 "latitude": -10.3333332,
    #                 "longitude": -53.1999999,
    #             }
    #         }
    #     }
    #     Note: geo reconstruction can cause disagreements between region and country attrs on internal nodes
    #     '''
    #     geo_lookup_json = {}
    #
    #
    #
    #     self.define_latitude_longitude()
    #     if "region" in geo_attributes:
    #         region_to_lat_long = {}
    #         regions = self.tree.get_attr_list("region")
    #         for region in regions:
    #             try:
    #                 region_to_lat_long[region] = self.location_to_lat_long[region]
    #             except:
    #                 print("REGION %s IS MISSING"%region)
    #         geo_lookup_json["region"] = region_to_lat_long
    #     if "country" in geo_attributes:
    #         country_to_lat_long = {}
    #         countries = self.tree.get_attr_list("country")
    #         for country in countries:
    #             country_to_lat_long[country] = self.location_to_lat_long[country]
    #             geo_lookup_json["country"] = country_to_lat_long
    #     if "division" in geo_attributes:
    #         division_to_lat_long = {}
    #         divisions = self.tree.get_attr_list("division")
    #         for division in divisions:
    #             division_to_lat_long[division] = self.location_to_lat_long[division]
    #         geo_lookup_json["division"] = division_to_lat_long
    #     return geo_lookup_json



    def auspice_export(self):
        '''
        export the tree, sequences, frequencies to json files for auspice visualization
        '''
        prefix = os.path.join(self.config["output"]["auspice"], self.info["prefix"])
        indent = 2
        # export json file that contains alignment diversity column by column
        self.seqs.export_diversity(fname=prefix+'_entropy.json', indent=2)
        # exports the tree and the sequences inferred for all clades in the tree
        if hasattr(self, 'tree') and self.tree is not None:
            self.tree.export(path=prefix, extra_attr = self.config["auspice"]["extra_attr"]
                         + ["muts", "aa_muts","attr", "clade"], indent = indent)


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
            write_json(freq_json, prefix+'_frequencies.json', indent=indent)

        # count number of tip nodes
        virus_count = 0
        for node in self.tree.tree.get_terminals():
            virus_count += 1

        # write out metadata json# Write out metadata
        print("Writing out metadata")
        meta_json = {}

        # join up config color options with those in the input JSONs.
        col_opts = self.config["auspice"]["color_options"]
        if self.colors:
            for trait, data in self.colors.iteritems():
                if trait in col_opts:
                    col_opts[trait]["color_map"] = [[k, v] for k, v in data.iteritems()]
                else:
                    self.log.warn("{} in colors (input JSON) but not auspice/color_options. Ignoring".format(trait))

        meta_json["color_options"] = col_opts
        meta_json["date_range"] = self.config["auspice"]["date_range"]
        meta_json["panels"] = self.config["auspice"]["panels"]
        meta_json["updated"] = time.strftime("X%d %b %Y").replace('X0','X').replace('X','')
        meta_json["virus_count"] = virus_count
        try:
            from pygit2 import Repository, discover_repository
            current_working_directory = os.getcwd()
            repository_path = discover_repository(current_working_directory)
            repo = Repository(repository_path)
            commit_id = repo[repo.head.target].id
            meta_json["commit"] = str(commit_id)
        except ImportError:
            meta_json["commit"] = "unknown"
        if len(self.config["auspice"]["controls"]):
            meta_json["controls"] = self.make_control_json(self.config["auspice"]["controls"])
        meta_json["geo"] = self.lat_longs
        write_json(meta_json, prefix+'_meta.json')

    def run_geo_inference(self):
        # run geo inference for all the things we have lat longs for
        if not self.lat_longs or len(self.lat_longs)==0:
            self.log.notify("no geo inference - no specified lat/longs")
            return
        for geo_attr in self.config["geo_inference"]:
            self.log.notify("running geo inference for {}".format(geo_attr))
            self.tree.geo_inference(geo_attr)

    def save_as_nexus(self):
        save_as_nexus(self.tree.tree, self.output_path + "_timeTree.nex")

if __name__=="__main__":
    print("This shouldn't be called as a script.")
