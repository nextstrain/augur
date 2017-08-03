from __future__ import division, print_function
import sys, os, time, gzip, glob
from base.io_util import make_dir, remove_dir, tree_to_json, write_json, myopen
import json
# from pdb import set_trace
from collections import defaultdict

def export_frequency_json(self, prefix, indent):

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
    else:
        self.log.notify("Cannot export frequencies - pivots do not exist")


def summise_publications_from_tree(tree):
    info = defaultdict(lambda: {"n": 0, "title": "?"})
    for clade in tree.find_clades():
        if "authors" in clade.attr:
            info[clade.attr["authors"]]["n"] += 1
            if "title" in clade.attr:
                info[clade.attr["authors"]]["title"] = clade.attr["title"]
    return info

def export_metadata_json(self, prefix, indent):
    self.log.notify("Writing out metadata")
    meta_json = {}

    # count number of tip nodes
    virus_count = 0
    for node in self.tree.tree.get_terminals():
        virus_count += 1
    meta_json["virus_count"] = virus_count


    meta_json["author_info"] = summise_publications_from_tree(self.tree.tree)

    # join up config color options with those in the input JSONs.
    col_opts = self.config["auspice"]["color_options"]
    if self.colors:
        for trait, data in self.colors.iteritems():
            if trait in col_opts:
                col_opts[trait]["color_map"] = data
            else:
                self.log.warn("{} in colors (input JSON) but not auspice/color_options. Ignoring".format(trait))

    meta_json["color_options"] = col_opts
    if "date_range" in self.config["auspice"]:
        meta_json["date_range"] = self.config["auspice"]["date_range"]
    if "analysisSlider" in self.config["auspice"]:
        meta_json["analysisSlider"] = self.config["auspice"]["analysisSlider"]
    meta_json["panels"] = self.config["auspice"]["panels"]
    meta_json["updated"] = time.strftime("X%d %b %Y").replace('X0','X').replace('X','')

    if "defaults" in self.config["auspice"]:
        meta_json["defaults"] = self.config["auspice"]["defaults"]

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
