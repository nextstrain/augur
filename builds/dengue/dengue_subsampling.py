from __future__ import division, print_function
import numpy as np
from collections import defaultdict
from pprint import pprint
from base.io_util import myopen
from base.titer_model import TiterCollection

def populate_counts(obj):
    sequence_count_total = defaultdict(int)
    sequence_count_region = defaultdict(int)
    seqs_to_count = obj.seqs.values()
    for seq in seqs_to_count:
        sequence_count_total[(seq.attributes['date'].year,
                              seq.attributes['date'].month)]+=1
        sequence_count_region[(seq.attributes['region'],
                              seq.attributes['date'].year,
                              seq.attributes['date'].month)]+=1
    return (sequence_count_total, sequence_count_region)

def dengue_subsampling(params, years_back, titer_values, force_include = []):
    ### Category: bin on region, year, month
    category = lambda x: (x.attributes['region'],
                          x.attributes['date'].year,
                          x.attributes['date'].month)

    ### Priority: # titer measurements, # unambiguous sites
    if titer_values is not None:
        titer_count = TiterCollection.count_strains(titer_values)

        def priority(seq):
            strain = seq.attributes['strain']
            accession = seq.attributes['accession']

            if strain in force_include or accession in force_include:
                return 10000

            if strain in titer_count:
                pr = titer_count[strain]
            else:
                pr = 0
            return pr + len(seq.seq)*0.00005 - 0.01*np.sum([seq.seq.count(nuc) for nuc in 'NRWYMKSHBVD'])

    else:
        print("Couldn't load titer information - using random priorities")
        def priority(seq):
            strain = seq.attributes['strain']
            accession = seq.attributes['accession']
            if strain in force_include or accession in force_include:
                return 10000
            else:
                return np.random.random()

    return {
        "category": category,
        "priority": priority,
        "threshold": 3
    }


    ### Per-region threshold
    # sampling_threshold = 3
    # region_threshold = int(np.ceil(1.0*sampling_threshold/len(regions)))
    # if type_of_subsampling == "even":
    #     def threshold(obj):
    #         """
    #         a higher order function which returns a fn which has access to
    #         some summary stats about the sequences (closure)
    #         """
    #         sequence_count_total, sequence_count_region = populate_counts(obj)
    #         def threshold_fn(x):
    #             #x is the collection key, in this case a tuple of (region, year, month)
    #             if sequence_count_total[(x[1], x[2])] < sampling_threshold:
    #                 return sampling_threshold
    #             region_counts = sorted([sequence_count_region[(r[0], x[1], x[2])] for r in regions])
    #             if region_counts[0] > region_threshold:
    #                 return region_threshold
    #             left_to_fill = sampling_threshold - len(regions)*region_counts[0]
    #             thres = region_counts[0]
    #             for ri, rc in zip(range(len(regions)-1, 0, -1), region_counts[1:]):
    #                 if left_to_fill - ri*(rc-thres)>0:
    #                     left_to_fill-=ri*(rc-thres)
    #                     thres = rc
    #                 else:
    #                     thres += left_to_fill/ri
    #                     break
    #             return max(1, int(thres))
    #         return threshold_fn
