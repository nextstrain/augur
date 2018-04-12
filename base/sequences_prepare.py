'''
parse, filter, subsample and save as JSON
'''
from __future__ import division, print_function
import os, re, time, csv, sys
from base.titer_model import TiterCollection
from io_util import myopen
from collections import defaultdict
from Bio import SeqIO
import numpy as np
from seq_util import pad_nucleotide_sequences, nuc_alpha, aa_alpha
from datetime import datetime
import json
from pdb import set_trace
import git
from utils import fix_names, num_date, ambiguous_date_to_date_range, potentially_combine
from pprint import pprint
import logging

logger = logging.getLogger(__name__)
TINY = 1e-10

class sequence_set(object):
    """ sequence set deals with loading sequences (stored in self.seqs)
    and various basic tasks including filtering, output etc """

    @classmethod
    def get_gene_name(cls, gene, genes):
        """Return a gene name for the given gene identifier and gene config.

        Args:
            cls (Python class): class for this method
            gene (str): gene identifier from a GenBank record
            genes (list or dict): a list of GenBank gene identifiers or a
                                  dictionary indexed by identified to a
                                  preferred gene name

        Returns:
            str: a gene name for the given gene identifier

        >>> genes = ["HA1", "HA2"]
        >>> sequence_set.get_gene_name("HA1", genes)
        'HA1'
        >>> genes = {"HA1": "HA1", "HA": "SigPep"}
        >>> sequence_set.get_gene_name("HA1", genes)
        'HA1'
        >>> sequence_set.get_gene_name("HA", genes)
        'SigPep'
        >>> sequence_set.get_gene_name("HA2", genes)
        'HA2'
        """
        try:
            return genes.get(gene, gene)
        except AttributeError:
            return gene

    def __init__(self, logger, segmentName):
        super(sequence_set, self).__init__()
        self.log = logger
        self.segmentName = segmentName
        self.extras = {}

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

    def _parse_date_per_seq(self, seq, fmts):
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
        # ## helpful debugging statements:
        # try:
        #     self.log.debug("{} date: {}, raw_date: {}, num_date: {}".format(seq.name, seq.attributes['date'], seq.attributes['raw_date'], seq.attributes['num_date']))
        # except KeyError:
        #     self.log.debug("{} missing date(s)".format(seq.name))

    def parse_date(self, fmts, prune):
        if not hasattr(self.seqs.values()[0], "attributes"):
            self.log.fatal("parse meta info first")
        from datetime import datetime
        for seq in self.seqs.values():
            self._parse_date_per_seq(seq, fmts)
        if prune:
            count = len(self.seqs)
            self.filterSeqs("Missing Date", lambda x:'date' in x.attributes and type(x.attributes['date'])!=str)
            self.log.notify("Removed sequences with missing dates (segment {}). n: {} -> {}".format(self.segmentName, count, len(self.seqs)))

    def filterSeqs(self, funcName, func):
        names = set(self.seqs.keys())
        if hasattr(self, "reference") and self.reference.include == 2:
            self.seqs = {key:seq for key, seq in self.seqs.iteritems() if func(seq) or key==self.reference.name}
        else:
            self.seqs = {key:seq for key, seq in self.seqs.iteritems() if func(seq)}
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
        # print("subsampling fn for ", name)
        try:
            fnIn = config["subsample"][name]
        except KeyError:
            # print("not found - using defaults")
            return defaults[name]

        # first - catch the case of INT threshold
        if name == "threshold" and type(fnIn) is int:
            # print("thresholding to INT", fnIn)
            return lambda x: fnIn

        # fnIn may be a higher order function
        try:
            fn = fnIn(self)
            assert(callable(fn))
            # print("HIGHER ORDER :)")
            return fn
        except Exception as e:
            pass
        try:
            assert(callable(fnIn))
            # print("CALLABLE :)")
            return fnIn
        except AssertionError:
            pass
        # print("not higher order or callable - falling back to defaults")
        return defaults[name]


    def get_subsampled_names(self, config):
        '''
        see docs/prepare.md for explination
        '''
        self.sequence_categories = defaultdict(list)
        names_prior = set(self.seqs.keys())
        seqs_to_subsample = self.seqs.values()

        if config.get("strains") is not None:
            with open(config["strains"], "r") as fh:
                include = set([line.strip() for line in fh])
        else:
            category = self.getSubsamplingFunc(config, "category")
            priority = self.getSubsamplingFunc(config, "priority")
            threshold = self.getSubsamplingFunc(config, "threshold")

            # sort sequences into categories and assign priority score
            for seq in seqs_to_subsample:
                seq._priority = priority(seq)
                self.sequence_categories[category(seq)].append(seq)

            # collect taxa to include
            include = set([])
            for cat, seqs in self.sequence_categories.iteritems():
                under_sampling = min(1.00, 1.0*len(seqs)/threshold(cat))
                for s in seqs: s.under_sampling=under_sampling
                seqs.sort(key=lambda x:x._priority, reverse=True)
                include.update([seq.id for seq in seqs[:threshold(cat)]])

        self.log.notify("Subsampling will take {}/{} strains".format(len(include), len(seqs_to_subsample)))
        return include

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
            "traits_are_dates": [],
            "title": config["title"],
            "maintainer": config["maintainer"],
            "auspice_filters": config["auspice_filters"]
        }
        if "traits_are_dates" in config and isinstance(config["traits_are_dates"], (list, tuple)):
            data["info"]["traits_are_dates"] = [trait for trait in config["traits_are_dates"] if trait in config["header_fields"].values()]
        data["info"]["prefix"] = prefix
        if self.segmentName == "genome":
            data["info"]["input_file"] = config["input_paths"][0]
        else:
            data["info"]["input_file"] = config["input_paths"][config["segments"].index(self.segmentName)]
        if "time_interval" in config:
            data["info"]["time_interval"] = [str(x) for x in config["time_interval"]]
        potentially_combine(config, data["info"], "regions")
        potentially_combine(config, data["info"], "lineage", False)
        data["info"]["segment"] = self.segmentName
        potentially_combine(config, data["info"], "resolution", False)
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
                # "attributes": self.reference.attributes,
                "strain": self.reference.attributes["strain"],
                "seq": str(self.reference.seq),
                "genes": self.reference.genes,
                "included": self.reference.name in self.seqs
            }

        # Titers must be present in the config and not None to be used.
        if config.get("titers") is not None:
            # Subset titer data to match the strains selected for export.
            filtered_titers = TiterCollection.filter_strains(config["titers"], self.seqs.keys())

            # Convert tuple dictionary keys to strings for JSON compatability.
            data["titers"] = {str(key): value
                              for key, value in filtered_titers.iteritems()}
            logger.debug("Filtered titers from %i to %i measures" % (len(config["titers"]), len(data["titers"])))

        # Flu-specific elements...
        if "vaccine_choices" in config and config["vaccine_choices"] is not None:
            data["info"]["vaccine_choices"] = {}
            for k, v in config["vaccine_choices"].items():
                if k in self.extras["leaves"]:
                    data["info"]["vaccine_choices"][k] = v
                else:
                    print("WARNING! Vaccine strain {} was not present in the data".format(k))
        if "LBI_params" in config:
            data["info"]["LBI_params"] = config["LBI_params"]
        if "frequency_params" in config:
            data["info"]["frequency_params"] = config["frequency_params"]

        json.dump(data, fh, indent=2)

    def load_reference(self, path, fmts, metadata, include=2, genes=False):
        """Assume it's genbank."""
        try:
            self.reference = SeqIO.read(path, 'genbank')
        except Exception as e:
            self.log.fatal("Problem reading reference {}. Error: {}".format(path, e))

        ## some checks
        try:
            assert("strain" in metadata)
            if include > 0:
                assert("date" in metadata)
        except AssertionError as e:
            self.log.fatal("Poorly defined reference. Error:".format(e))

        if genes:
            # we used to make these FeatureLocation objects here, but that won't go to JSON
            # so just do it in the Process part instead. For reference:
            # FeatureLocation(start=f.location.start, end=f.location.end, strand=1)
            self.reference.genes = {
                sequence_set.get_gene_name(f.qualifiers['gene'][0], genes): {"start": int(f.location.start),
                                           "end": int(f.location.end), "strand": f.location.strand}
                for f in self.reference.features
                if 'gene' in f.qualifiers and f.qualifiers['gene'][0] in genes
            }
        else:
            self.reference.genes = {}

        # use the supplied metadata dict to define attributes
        seq_attr_keys = self.seqs.values()[0].attributes.keys()
        self.reference.attributes = {k:fix_names(v) for k,v in metadata.items() if k in seq_attr_keys}
        self.reference.name = self.reference.attributes["strain"]
        self.reference.id = self.reference.attributes["strain"]

        # is there any possibility that the reference will be added to the sequences?
        self.reference.include = include; # flag {0,1,2}
        if self.reference.name in self.seqs:
            self.log.notify("Segment {} reference already in dataset".format(self.segmentName))
            if include == 0:
                self.log.notify("Removing reference from pool of sequences to analyse")
                del self.seqs[self.reference.name]
        elif include > 0:
            ## add to sequences (tidy up attributes first)
            self._parse_date_per_seq(self.reference, fmts)
            self.seqs[self.reference.name] = self.reference
            missing_attrs = set(seq_attr_keys) - set(self.reference.attributes.keys()) - set(["date", "num_date"])
            if len(missing_attrs) > 0:
                self.log.notify("Including reference in segment {} but the following attributes are missing: {}".format(self.segmentName, " & ".join(missing_attrs)))

if __name__=="__main__":
    pass
