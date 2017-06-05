'''
parse, filter, subsample and save as JSON
'''
from __future__ import division, print_function
import os, re, time, csv, sys
from io_util import myopen
from collections import defaultdict
from Bio import SeqIO
import numpy as np
from seq_util import pad_nucleotide_sequences, nuc_alpha, aa_alpha
from datetime import datetime
import json
from pdb import set_trace
import git
from utils import fix_names, num_date, ambiguous_date_to_date_range
from pprint import pprint

TINY = 1e-10

class sequence_set(object):
    """ sequence set deals with loading sequences (stored in self.seqs)
    and various basic tasks including filtering, output etc """

    def __init__(self, logger, segmentName):
        super(sequence_set, self).__init__()
        self.log = logger
        self.segmentName = segmentName
        self.extras = {}
        self.reference = False

    def load_mfa(self, path):
        try:
            with myopen(path) as seq_file:
                self.seqs = {x.name:x for x in SeqIO.parse(seq_file, 'fasta')}
        except Exception as e:
            self.log.fatal("Error loading sequences from {}. Error: {}".format(path, e))
        self.nstart = len(self.seqs)
        self.log.notify("Loaded {} sequences from {}".format(self.nstart, path))

    def ungap(self):
        '''
        remove previously existing gaps and make sure all sequences are upper case
        '''
        for seq in self.seqs.values():
            seq.seq = seq.seq.ungap('-').upper()

    def parse_headers(self, fields, sep='|', strip='_'):
        '''
        split the sequence description and add annotations to sequences
        '''
        try:
            assert("strain" in fields.values())
        except AssertionError:
            self.log.fatal("Config file: fasta_fields must contain 'strain'")
        for seq in self.seqs.values():
            if not hasattr(seq, "attributes"): seq.attributes = {}
            words = map(lambda x:fix_names(x), seq.description.replace(">","").split(sep))
            for ii, val in enumerate(words):
                if ii in fields:
                    if val not in ["", "-"]:
                        # self.log.debug("{} -> {}".format(fields[ii], val))
                        seq.attributes[fields[ii]] = val
                    else:
                        seq.attributes[fields[ii]] = ""
        self.seqs = {seq.attributes['strain']:seq for seq in self.seqs.values()}
        for seq in self.seqs.values():
            seq.id = seq.attributes['strain']
            seq.name = seq.attributes['strain']

    def load_json(self, path):
        pass

    def parse_date(self, fmts, prune):
        if not hasattr(self.seqs.values()[0], "attributes"):
            self.log.fatal("parse meta info first")
        from datetime import datetime
        for seq in self.seqs.values():
            if 'date' in seq.attributes and seq.attributes['date']!='':
                for fmt in fmts:
                    try:
                        if 'XX' in seq.attributes['date']:
                            min_date, max_date = ambiguous_date_to_date_range(seq.attributes['date'], fmt)
                            seq.attributes['raw_date'] = seq.attributes['date']
                            seq.attributes['num_date'] = np.array((num_date(min_date), num_date(max_date)))
                            seq.attributes['date'] = min_date
                        else:
                            if callable(fmt):
                                tmp = fmt(seq.attributes['date'])
                            else:
                                try:
                                    tmp = datetime.strptime(seq.attributes['date'], fmt).date()
                                except:
                                    tmp = seq.attributes['date']
                            seq.attributes['raw_date'] = seq.attributes['date']
                            seq.attributes['num_date'] = num_date(tmp)
                            seq.attributes['date']=tmp
                            break
                    except:
                        continue

        # helpful debugging statements:
        # for seq in self.seqs.values():
        #     try:
        #         self.log.debug("{} date: {}, raw_date: {}, num_date: {}".format(seq.name, seq.attributes['date'], seq.attributes['raw_date'], seq.attributes['num_date']))
        #     except KeyError:
        #         self.log.debug("{} missing date(s)".format(seq.name))

        if prune:
            self.filterSeqs("Missing Date", lambda x:'date' in x.attributes and type(x.attributes['date'])!=str)

    def filterSeqs(self, funcName, func):
        names = set(self.seqs.keys())
        self.seqs = {key:seq for key, seq in self.seqs.iteritems() if func(seq)} #or key==self.reference_seq.name
        for name in names - set(self.seqs.keys()):
            self.log.drop(name, self.segmentName, funcName)

    # def filter(self, func, leave_ref=False):
    #     if leave_ref:
    #         self.all_seqs = {key:seq for key, seq in self.all_seqs.iteritems() if func(seq) or key==self.reference_seq.name}
    #     else:
    #         self.all_seqs = {key:seq for key, seq in self.all_seqs.iteritems() if func(seq)}
    #     print("Filtered to %d sequences"%len(self.all_seqs))

    def getSubsamplingFunc(self, config, name):
        defaults = {
            "category": lambda x: (x.attributes['date'].year, x.attributes['date'].month),
            "priority": lambda x: np.random.random(),
            "threshold": lambda x: 5
        }
        print("subsampling fn for ", name)
        try:
            fnIn = config["subsample"][name]
        except KeyError:
            print("not found - using defaults")
            return defaults[name]

        # first - catch the case of INT threshold
        if name == "threshold" and type(fnIn) is int:
            print("thresholding to INT", fnIn)
            return lambda x: fnIn

        # fnIn may be a higher order function
        try:
            fn = fnIn(self)
            assert(callable(fn))
            print("HIGHER ORDER :)")
            return fn
        except Exception as e:
            pass
        try:
            assert(callable(fnIn))
            print("CALLABLE :)")
            return fnIn
        except AssertionError:
            pass
        print("not higher order or callable - falling back to defaults")
        return defaults[name]


    def subsample(self, config):
        '''
        see docs/prepare.md for explination
        '''
        category = self.getSubsamplingFunc(config, "category")
        priority = self.getSubsamplingFunc(config, "priority")
        threshold = self.getSubsamplingFunc(config, "threshold")

        self.sequence_categories = defaultdict(list)
        names_prior = set(self.seqs.keys())
        seqs_to_subsample = self.seqs.values()

        # sort sequences into categories and assign priority score
        for seq in seqs_to_subsample:
            seq._priority = priority(seq)
            self.sequence_categories[category(seq)].append(seq)

        # sample and record the degree to which a category is under_sampled
        self.seqs = {}
        for cat, seqs in self.sequence_categories.iteritems():
            under_sampling = min(1.00, 1.0*len(seqs)/threshold(cat))
            for s in seqs: s.under_sampling=under_sampling
            seqs.sort(key=lambda x:x._priority, reverse=True)
            self.seqs.update({seq.id:seq for seq in seqs[:threshold(cat)]})

        self.log.notify("Subsampling segment {}. n={} -> {}".format(self.segmentName, len(seqs_to_subsample), len(self.seqs)))
        for name in names_prior - set(self.seqs.keys()):
            self.log.drop(name, self.segmentName, "subsampled")


    # def strip_non_reference(self):
    #     ungapped = np.array(self.sequence_lookup[self.reference_seq.name])!='-'
    #     from Bio.Seq import Seq
    #     for seq in self.aln:
    #         seq.seq = Seq("".join(np.array(seq)[ungapped]))

    def get_trait_values(self, trait):
        vals = set()
        for seq, obj in self.seqs.iteritems():
            if trait in obj.attributes:
                vals.add(obj.attributes[trait])
        # don't forget the reference here
        try:
            if self.reference.use:
                vals.add(self.reference.attributes[trait])
        except:
            pass # perhaps not set, no refrence, whatever
        return vals

    def write_json(self, fh, config, prefix):
        # datetime() objects and [arrays] don't go to JSONs
        # not a problem - we still have raw_date to get them back
        for seq in self.seqs.values():
            if 'date' in seq.attributes:
                del seq.attributes['date']
            if 'num_date' in seq.attributes:
                del seq.attributes['num_date']

        data = self.extras
        data["info"] = {
            "n(starting)": self.nstart,
            "n(final)": len(self.seqs),
            "commit": git.Repo(search_parent_directories=True).head.object.hexsha,
            "date_format": config["date_format"],
            "subsampled": bool(config["subsample"]),
            "traits_are_dates": []
        }
        if "traits_are_dates" in config and isinstance(config["traits_are_dates"], (list, tuple)):
            data["info"]["traits_are_dates"] = [trait for trait in config["traits_are_dates"] if trait in config["header_fields"].values()]
        data["info"]["prefix"] = prefix
        if self.segmentName == "genome":
            data["info"]["input_file"] = config["input_paths"][0]
        else:
            data["info"]["input_file"] = config["input_paths"][config["segments"].index(self.segmentName)]
        if config["time_interval"]:
            data["info"]["time_interval"] = [str(x) for x in config["time_interval"]]
        if config["regions"]:
            data["info"]["regions"] = config["regions"]


        data["sequences"] = {}
        for seqName, seq in self.seqs.iteritems():
            data["sequences"][seqName] = {
                "attributes": seq.attributes,
                "seq": str(seq.seq)
            }

        if not self.reference:
            data["reference"] = None
        else:
            data["reference"] = {
                "attributes": self.reference.attributes,
                "seq": str(self.reference.seq),
                "genes": self.reference.genes
            }


        json.dump(data, fh, indent=2)

    def load_reference(self, path, metadata, use=True, genes=False):
        """Assume it's genbank."""
        try:
            self.reference = SeqIO.read(path, 'genbank')
        except Exception as e:
            self.log.fatal("Problem reading reference {}. Error: {}".format(path, e))

        ## some checks
        try:
            assert("strain" in metadata)
        except AssertionError as e:
            self.log.fatal("Poorly defined reference. Error:".format(e))

        # use the metadata to define things!
        seq_attr_keys = self.seqs.values()[0].attributes.keys()
        self.reference.attributes = {k:fix_names(v) for k,v in metadata.items() if k in seq_attr_keys}
        # date gets changed to raw_date for sequences. be consistent.
        if 'date' in self.reference.attributes.keys():
            self.reference.attributes["raw_date"] = self.reference.attributes["date"]
            del self.reference.attributes["date"]

        if genes:
            # we used to make these FeatureLocation objects here, but that won't go to JSON
            # so just do it in the Process part instead. For reference:
            # FeatureLocation(start=f.location.start, end=f.location.end, strand=1)
            self.reference.genes = {f.qualifiers['gene'][0]:{"start": int(f.location.start), "end": int(f.location.end), "strand": 1} for f in self.reference.features if 'gene' in f.qualifiers and f.qualifiers['gene'][0] in genes}
        else:
            self.reference.genes = {}

if __name__=="__main__":
    pass
