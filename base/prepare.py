from __future__ import division, print_function
import sys, os
import numpy as np
from datetime import datetime
from base.sequences_new import sequence_set
from base.logger import logger
from base.utils import generate_cmap, define_latitude_longitude
# useful for debugging
from pdb import set_trace
from pprint import pprint

required_config_fields = [
    "dir", "file_prefix", "segments", "input_format", "input_paths",
    "output_folder", "header_fields", "date_format", "ensure_all_segments",
    "require_dates", "colors", "lat_longs"
]

class prepare(object):
    def __init__(self, config):
        """ check config file, make necessary directories, set up logger """
        super(prepare, self).__init__()

        try:
            for x in required_config_fields:
                assert(x in config)
        except AssertionError:
            print("Fatal Error: Config file is missing field '{}'".format(x))
            sys.exit(2)

        try:
            assert(os.path.basename(os.getcwd()) == config["dir"])
        except AssertionError:
            print("Run this script from within the {} directory".format(config["dir"]))
            sys.exit(2)

        for p in [config["output_folder"]]:
            if not os.path.isdir(p):
                os.makedirs(p)

        self.log = logger(config["output_folder"], False)
        self.config = config

        try:
            if config["segments"] == False:
                assert(len(config["input_paths"]) == 1)
                config["segments"] = ["genome"]
            else:
                assert(len(config["segments"]) == len(config["input_paths"]))
        except AssertionError:
            self.log.fatal("Config file: # segments don't match # input paths")
        try:
            for p in config["input_paths"]:
                assert(os.path.exists(p))
        except AssertionError:
            self.log.fatal("Config file: input path '{}' doesn't exist".format(p))

        # this block initialses this.segments
        if config["input_format"] == "fasta":
            self.load_fasta()
        else:
            self.log.fatal("Currently only FASTA sequences can be loaded".format(config["dir"]))


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
        """ subsample only the first one, then filter the rest such that the names match """
        if "subsample" not in self.config or self.config["subsample"] == False:
            return
        # subsample the first (may be the only one)
        self.segments[self.config["segments"][0]].subsample(self.config)
        subsampled_taxa = self.segments[self.config["segments"][0]].seqs.keys()
        for idx in range(1, len(self.config["segments"])):
            self.segments[self.config["segments"][idx]].filterSeqs("subsampled", lambda s: s.id in subsampled_taxa)
            self.log.notify("Matched segment {} to subsampled taxa. n: {}".format(self.config["segments"][idx], len(self.segments[self.config["segments"][idx]].seqs)))


    def ensure_all_segments(self):
        if "ensure_all_segments" in self.config and self.config["ensure_all_segments"] == True:
            in_all = set.intersection(*[set(obj.seqs.keys()) for seg, obj in self.segments.iteritems()])
            # to_drop = {n for segment in self.segments.keys() for n in self.segments[segment].seqs.keys()} - in_all
            self.log.notify("Removing sequences without the full complement of segments")
            for seg, obj in self.segments.iteritems():
                obj.filterSeqs("ensure_all_segments", lambda s: s.name in  in_all)
            self.log.notify("n(seqs): "+ ", ".join([a+"="+str(len(b.seqs)) for a,b in self.segments.items()]))

    """for the proper english speakers"""
    def colours(self):
        self.colors()

    def colors(self):
        if self.config["colors"]:
            cols = {}
            for trait in self.config["colors"]:
                try:
                    assert(trait in self.config["header_fields"].values())
                except AssertionError:
                    self.log.warn("Colour trait {} not in header_fields. Skipping".format(trait))
                    continue
                self.log.notify("Generating colour maps for '{}'".format(trait))
                vals = set.union(*[obj.get_trait_values(trait) for seg, obj in self.segments.iteritems()])
                cols[trait] = generate_cmap(vals, False)

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
            for trait in self.config["colors"]:
                vals = set.union(*[obj.get_trait_values(trait) for seg, obj in self.segments.iteritems()])
                lat_longs[trait] = {}
                # do it the long way to handle missing data
                for key in vals: # well named
                    try:
                        lat_longs[trait][key] = lat_long_db[key]
                    except KeyError:
                        lat_longs[trait][key] = "unknown"
                        self.log.warn("Unknown lat/longs for {} {}. Attempting to continue but this will most likely cause problems in Auspice.".format(trait, key))

            # save to each sequence_set object. It's them that write the JSONs
            for seg, obj in self.segments.iteritems():
                obj.extras["lat_longs"] = lat_longs


    def write_to_json(self):
        for seg, obj in self.segments.iteritems():
            if seg == "genome":
                fname = os.path.join(self.config["output_folder"], self.config["file_prefix"]+".json")
            else:
                fname = os.path.join(self.config["output_folder"], self.config["file_prefix"]+"_"+seg+".json")
            with open(fname, 'w') as fh:
                obj.write_json(fh, self.config)
