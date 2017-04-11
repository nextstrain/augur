from __future__ import division, print_function
import sys, os
import numpy as np
from datetime import datetime
from base.sequences_new import sequence_set
from base.logger import logger
from pdb import set_trace

required_config_fields = [
    "dir", "file_prefix", "segments", "input_format", "input_paths",
    "output_folder", "header_fields", "date_format", "ensure_all_segments",
    "require_dates"
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
                    self.log.notify("Applying filter '{}' to segment '{}'".format(fName, seg))
                    obj.filterSeqs(fName, fFunc)
            else: # wasn't callable
                if isinstance(fFunc, dict):
                    for seg, segFunc in fFunc.iteritems():
                        if callable(segFunc) and seg in self.segments.keys():
                            self.log.notify("Applying filter '{}' to segment '{}'".format(fName, seg))
                            self.segments[seg].filterSeqs(fName, segFunc)
                        else:
                            self.log.warn("Filter {} not applied to segment {} (not callable or segment not found)".format(fName, seg))
                else:
                    self.log.warn("Filter {} not applied - neither callable or dict of callables".format(fName))

    def subsample(self):
        if "subsample" not in self.config or self.config["subsample"] == False:
            return
        for seg, obj in self.segments.iteritems():
            obj.subsample(self.config)

    def ensure_all_segments(self):
        if "ensure_all_segments" in self.config and self.config["ensure_all_segments"] == True:
            in_all = set.intersection(*[set(obj.seqs.keys()) for seg, obj in self.segments.iteritems()])
            # to_drop = {n for segment in self.segments.keys() for n in self.segments[segment].seqs.keys()} - in_all
            self.log.notify("Removing sequences without the full complement of segments")
            for seg, obj in self.segments.iteritems():
                obj.filterSeqs("ensure_all_segments", lambda s: s.name in  in_all)
            self.log.notify("n(seqs): "+ ", ".join([a+"="+str(len(b.seqs)) for a,b in self.segments.items()]))


    def write_to_json(self):
        for seg, obj in self.segments.iteritems():
            if seg == "genome":
                fname = os.path.join(self.config["output_folder"], self.config["file_prefix"]+".json")
            else:
                fname = os.path.join(self.config["output_folder"], self.config["file_prefix"]+"_"+seg+".json")
            with open(fname, 'w') as fh:
                obj.write_json(fh, self.config)
