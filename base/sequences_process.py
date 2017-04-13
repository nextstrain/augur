from __future__ import division, print_function
import os, re, time, csv, sys
# from io_util import myopen
# from io_util import myopen, make_dir, remove_dir, tree_to_json, write_json
# from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import FeatureLocation
import numpy as np
# from seq_util import pad_nucleotide_sequences, nuc_alpha, aa_alpha
from datetime import datetime
import json
from pdb import set_trace
from utils import fix_names, num_date, ambiguous_date_to_date_range
from pprint import pprint

class sequence_set(object):

    def __init__(self, logger, sequences, reference, dateFormat):
        super(sequence_set, self).__init__()
        self.log = logger
        self.reference = None
        self.proteins = None

        # load sequences from the (parsed) JSON - don't forget to sort out dates
        self.seqs = {}
        for name, data in sequences.iteritems():
            self.seqs[name] = SeqRecord(Seq(data["seq"], generic_dna),
                   id=name, name=name, description=name)
            self.seqs[name].attributes = data["attributes"]
            # tidy up dates
            self.parse_date(self.seqs[name], dateFormat)

        # load reference from (parsed) JSON & clean up dates
        if reference and len(reference):
            name = reference["attributes"]["strain"]
            self.reference_seq = SeqRecord(Seq(reference["seq"], generic_dna),
                   id=name, name=name, description=name)
            self.reference_seq.attributes = reference["attributes"]
            self.parse_date(self.reference_seq, dateFormat)
            # is reference already in self.seqs?

            #sort out the proteins:
            if "genes" in reference and len(reference["genes"]):
                self.proteins = {k:FeatureLocation(start=v["start"], end=v["end"], strand=v["strand"]) for k, v in reference["genes"].iteritems()}






    """ this function is similar to that in sequences_prepare
    these should be consolidated
    """
    def parse_date(self, seq, fmts):
        for fmt in fmts:
            try:
                if 'XX' in seq.attributes['raw_date']:
                    min_date, max_date = ambiguous_date_to_date_range(seq.attributes['raw_date'], fmt)
                    # seq.attributes['raw_date'] = seq.attributes['date']
                    seq.attributes['num_date'] = np.array((num_date(min_date), num_date(max_date)))
                    seq.attributes['date'] = min_date
                else:
                    if callable(fmt):
                        tmp = fmt(seq.attributes['raw_date'])
                    else:
                        try:
                            tmp = datetime.strptime(seq.attributes['raw_date'], fmt).date()
                        except:
                            tmp = seq.attributes['raw_date']
                    # seq.attributes['raw_date'] = seq.attributes['date']
                    seq.attributes['num_date'] = num_date(tmp)
                    seq.attributes['date']=tmp
                    break
            except:
                continue
