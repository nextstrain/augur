from __future__ import division, print_function
import argparse
import sys, os
import numpy as np
from datetime import datetime
from base.config import combine_configs
from base.sequences_prepare import sequence_set
from base.logger import logger
from base.utils import generate_cmap, define_latitude_longitude
# useful for debugging
from pdb import set_trace
from pprint import pprint
from collections import defaultdict

import logging
logging.basicConfig(
    format='[%(levelname)s] %(asctime)s - %(name)s - %(message)s',
    level=logging.INFO,
    stream=sys.stderr
)


def collect_args():
    """Returns a minimal default argument parser instance for all prepare scripts.
    """
    parser = argparse.ArgumentParser(
        description="Prepare fauna FASTA for analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-v', '--viruses_per_month', type = int, help='Subsample x viruses per country per month. Set to 0 to disable subsampling.')
    parser.add_argument('--sequences', nargs='+', help="FASTA file of virus sequences from fauna (e.g., zika.fasta)")
    parser.add_argument('--file_prefix', help="Prefix for output files (e.g., 'zika' for output like 'auspice/zika_meta.json')")
    parser.add_argument('--verbose', action="store_true", help="turn on verbose reporting")

    parser.set_defaults(viruses_per_month=15)

    return parser

class prepare(object):
    def __init__(self, config):
        """ check config file, make necessary directories, set up logger """
        super(prepare, self).__init__()
        self.config = combine_configs("prepare", config)

        try:
            assert(os.path.basename(os.getcwd()) == self.config["dir"])
        except AssertionError:
            print("Run this script from within the {} directory".format(self.config["dir"]))
            sys.exit(2)

        for p in [self.config["output_folder"]]:
            if not os.path.isdir(p):
                os.makedirs(p)

        self.log = logger(self.config["output_folder"], False)

        try:
            if self.config["segments"] == False:
                assert(len(self.config["input_paths"]) == 1)
                self.config["segments"] = ["genome"]
            else:
                assert(len(self.config["segments"]) == len(self.config["input_paths"]))
        except AssertionError:
            self.log.fatal("Config file: # segments don't match # input paths")
        try:
            for p in self.config["input_paths"]:
                assert(os.path.exists(p))
        except AssertionError:
            self.log.fatal("Config file: input path '{}' doesn't exist".format(p))

        # this block initialses this.segments
        if self.config["input_format"] == "fasta":
            self.load_fasta()
        else:
            self.log.fatal("Currently only FASTA sequences can be loaded".format(self.config["dir"]))

    def load_fasta(self):
        self.segments = {}
        for idx in range(0, len(self.config["segments"])):
            segmentName = self.config["segments"][idx]
            segmentObj = sequence_set(self.log, segmentName)
            segmentObj.load_mfa(self.config["input_paths"][idx])
            segmentObj.ungap()
            segmentObj.parse_headers(self.config["header_fields"])
            segmentObj.parse_date(self.config["date_format"], prune=self.config["require_dates"])
            self.segments[segmentName] = segmentObj

    def applyFilters(self):
        if self.config.get("strains") is not None:
            self.log.notify("Specific strains requested, skipping filters")
            return

        for (fName, fFunc) in self.config["filters"]:
            if callable(fFunc):
                for seg, obj in self.segments.iteritems():
                    tmp = len(obj.seqs)
                    obj.filterSeqs(fName, fFunc)
                    self.log.notify("Applied filter '{}' to segment '{}'. n: {} -> {}".format(fName, seg, tmp, len(obj.seqs)))
            else: # wasn't callable
                if isinstance(fFunc, dict):
                    for seg, segFunc in fFunc.iteritems():
                        if callable(segFunc) and seg in self.segments.keys():
                            tmp = len(self.segments[seg].seqs)
                            self.segments[seg].filterSeqs(fName, segFunc)
                            self.log.notify("Applied filter '{}' to segment '{}'. n: {} -> {}".format(fName, seg, tmp, len(self.segments[seg].seqs)))
                        else:
                            self.log.warn("Filter {} not applied to segment {} (not callable or segment not found)".format(fName, seg))
                else:
                    self.log.warn("Filter {} not applied - neither callable or dict of callables".format(fName))

    def subsample(self):
        """ subsample only the first segment, then filter the rest such that the names match
        note that if reference.include == 2 (per segment) there may be different seqs per segment
        """
        if "subsample" not in self.config or self.config["subsample"] == False:
            return
        taxa_to_include = self.segments[self.config["segments"][0]].get_subsampled_names(self.config)
        for name in self.config["segments"]:
            self.segments[name].filterSeqs("subsampled", lambda s: s.id in taxa_to_include)
            refStr = "(incl. ref)" if self.segments[name].reference.name in self.segments[name].seqs.keys() else "(no Ref)"
            self.log.notify("Subsampled segment {}. n={} {}".format(name, len(self.segments[name].seqs), refStr))

    def ensure_all_segments(self):
        if "ensure_all_segments" in self.config and self.config["ensure_all_segments"] == True:
            if len(self.segments) == 1:
                self.log.notify("Ignoring \"ensure_all_segments\" as there's only 1 segment...")
                return
            in_all = set.intersection(*[set(obj.seqs.keys()) for seg, obj in self.segments.iteritems()])
            # to_drop = {n for segment in self.segments.keys() for n in self.segments[segment].seqs.keys()} - in_all
            self.log.notify("Removing sequences without the full complement of segments")
            for seg, obj in self.segments.iteritems():
                obj.filterSeqs("ensure_all_segments", lambda s: s.name in  in_all)
                # NB if reference.include == 2 the refs will be kept (per segment)
            self.log.notify("n(seqs): "+ ", ".join([a+"="+str(len(b.seqs)) for a,b in self.segments.items()]))

    """for the proper english speakers"""
    def colours(self):
        self.colors()

    def parse_colors_file(self):
        db = defaultdict(list)
        if "color_defs" not in self.config: return db
        try:
            for fname in self.config["color_defs"]:
                with open(fname) as fh:
                    for line in fh:
                        # line: trait   trait_value     hex_code
                        if line.startswith('#'): continue
                        fields = line.strip().split()
                        if len(fields)!=3: continue
                        db[fields[0]].append((fields[1], fields[2]))
        except IOError:
            self.log.warn("Couldn't open color definitions file {}.".format(self.config["color_defs"]))
        return db

    def colors(self):
        if self.config["colors"]:
            cols = self.parse_colors_file() #defaultdict type <list>
            for trait in self.config["colors"]:
                try:
                    assert(trait in self.config["header_fields"].values())
                except AssertionError:
                    self.log.warn("Colour trait {} not in header_fields. Skipping".format(trait))
                    continue
                self.log.notify("Generating colour maps for '{}'".format(trait))
                vals = set.union(*[obj.get_trait_values(trait) for seg, obj in self.segments.iteritems()])
                missing_vals = vals - set([x[0] for x in cols[trait]]) - set(["?"])
                if len(missing_vals):
                    cols[trait].extend(generate_cmap(missing_vals, False))
                # remove colours for values that have no sequences
                for idx in [i for i, v in enumerate(cols[trait]) if v[0] not in vals][::-1]:
                    del cols[trait][idx]
            # save to each sequence_set object. It's them that write the JSONs
            for seg, obj in self.segments.iteritems():
                obj.extras["colors"] = cols

    def latlongs(self):
        if self.config["lat_longs"]:
            lat_longs = {}
            try:
                assert("lat_long_defs" in self.config)
            except AssertionError:
                self.log.fatal("You asked for lat/longs but didn't provide definition files ('lat_long_defs')")
            lat_long_db = define_latitude_longitude(self.config["lat_long_defs"], self.log)
            if not len(lat_long_db):
                self.log.fatal("You asked for lat/longs but the definition files weren't useful")
            for trait in self.config["lat_longs"]:
                vals = set.union(*[obj.get_trait_values(trait) for seg, obj in self.segments.iteritems()])
                lat_longs[trait] = {}
                # do it the long way to handle missing data
                for key in vals: # well named
                    try:
                        lat_longs[trait][key] = lat_long_db[key.lower()]
                    except KeyError:
                        if key != "?":
                            lat_longs[trait][key] = {'latitude': 0,'longitude': 0}
                            self.log.warn("Unknown lat/longs for {} {}. Setting to 0,0 in order to appease auspice but you should fix this.".format(trait, key))

            # save to each sequence_set object. It's them that write the JSONs
            for seg, obj in self.segments.iteritems():
                obj.extras["lat_longs"] = lat_longs

    def load_references(self):
        if "reference" in self.config:
            if len(self.config["segments"]) > 1:
                self.log.fatal("Config must define references (not reference) for segemented virus")
            self.segments[self.config["segments"][0]].load_reference(fmts=self.config["date_format"], **self.config["reference"])
        elif "references" in self.config:
            for k, v in self.config["references"].iteritems():
                if k not in self.segments.keys():
                    self.log.warn("Reference segment {} not found in segments".format(k))
                else:
                    self.segments[k].load_reference(fmts=self.config["date_format"], **v)

    def write_to_json(self, segment_addendum=False):
        for seg, obj in self.segments.iteritems():
            prefix = self.config["file_prefix"]
            if segment_addendum:                                # add segment to file
                if "*segment*" in prefix:                       # replace *segment* placeholder if it exists
                    prefix = prefix.replace("*segment*", seg)
                else:                                           # else, append
                    prefix += "_" + seg
            fname = os.path.join(self.config["output_folder"], prefix + ".json")
            with open(fname, 'w') as fh:
                obj.write_json(fh, self.config, prefix)
