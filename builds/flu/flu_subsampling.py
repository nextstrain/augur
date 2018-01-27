from __future__ import division, print_function
from flu_info import regions,reference_viruses
import numpy as np
from collections import defaultdict
from pprint import pprint
from base.io_util import myopen
from base.titer_model import TiterCollection

vpm_dict = {
    2: 92,
    3: 62,
    6: 32,
    12: 18,
}

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

def flu_subsampling(params, years_back, titer_values):
    if params.sampling == "even":
        type_of_subsampling = "even"
    elif params.sampling in [x[0] for x in regions]:
        type_of_subsampling = "priority"
    else:
        type_of_subsampling = "flat"

    #### DEFINE THE CATEGORY:
    if type_of_subsampling in ["even", "priority"]:
        category = lambda x: (x.attributes['region'],
                              x.attributes['date'].year,
                              x.attributes['date'].month)
    else:
        category = lambda x: (x.attributes['date'].year,
                              x.attributes['date'].month)

    #### DEFINE THE PRIORITY
    if titer_values is not None:
        HI_titer_count = TiterCollection.count_strains(titer_values)
    else:
        print("Couldn't load titer information - using random priorities")
        HI_titer_count = False
        def priority(seq):
            return np.random.random() + int(seq.name in reference_viruses[params.lineage])
    if HI_titer_count:
        def priority(seq):
            sname = seq.attributes['strain']
            if sname in HI_titer_count:
                pr = HI_titer_count[sname]
            else:
                pr = 0
            return (pr + len(seq.seq)*0.0001 - 0.01*np.sum([seq.seq.count(nuc) for nuc in 'NRWYMKSHBVD']) +
                    1e6*int(seq.name in reference_viruses[params.lineage]))

    ##### DEFINE THE THRESHOLD
    if params.viruses_per_month != 0:
        sampling_threshold = params.viruses_per_month
    else:
        sampling_threshold = vpm_dict[years_back]

    region_threshold = int(np.ceil(1.0*sampling_threshold/len(regions)))
    if type_of_subsampling == "priority":
        priority_region = params.sampling
    if type_of_subsampling == "even":
        def threshold(obj):
            """
            a higher order function which returns a fn which has access to
            some summary stats about the sequences (closure)
            """
            sequence_count_total, sequence_count_region = populate_counts(obj)
            def threshold_fn(x):
                #x is the collection key, in this case a tuple of (region, year, month)
                if sequence_count_total[(x[1], x[2])] < sampling_threshold:
                    return sampling_threshold
                region_counts = sorted([sequence_count_region[(r[0], x[1], x[2])] for r in regions])
                if region_counts[0] > region_threshold:
                    return region_threshold
                left_to_fill = sampling_threshold - len(regions)*region_counts[0]
                thres = region_counts[0]
                for ri, rc in zip(range(len(regions)-1, 0, -1), region_counts[1:]):
                    if left_to_fill - ri*(rc-thres)>0:
                        left_to_fill-=ri*(rc-thres)
                        thres = rc
                    else:
                        thres += left_to_fill/ri
                        break
                return max(1, int(thres))
            return threshold_fn
    elif type_of_subsampling == "priority":
        priority_region = params.sampling
        fraction = 0.5
        def threshold(obj):
            """
            a higher order function which returns a fn which has access to
            some summary stats about the sequences (closure)
            """
            sequence_count_total, sequence_count_region = populate_counts(obj)
            def threshold_fn(x):
                #x is the collection key, in this case a tuple of (region, year, month)
                if x[0]==priority_region:
                    return int(sampling_threshold*fraction)
                nregions = len(regions)-1
                total_threshold_world = sampling_threshold*(1-fraction)
                region_threshold = int(np.ceil(1.0*total_threshold_world/nregions))
                region_counts = sorted([sequence_count_region[(r[0], x[1], x[2])]
                                        for r in regions if r!=priority_region])
                if region_counts[0]>region_threshold:
                    return region_threshold
                else:
                    left_to_fill = total_threshold_world - nregions*region_counts[0]
                    thres = region_counts[0]
                    for ri, rc in zip(range(nregions-1, 0, -1), region_counts[1:]):
                        if left_to_fill - ri*(rc-thres)>0:
                            left_to_fill-=ri*(rc-thres)
                            thres = rc
                        else:
                            thres += left_to_fill/ri
                            break
                    return max(1,int(thres))
            return threshold_fn
    else: # flat subsampling
        threshold = lambda x: sampling_threshold

    return {
        "category": category,
        "priority": priority,
        "threshold": threshold
    }
