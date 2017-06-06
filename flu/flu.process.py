from __future__ import print_function
import os, sys, glob
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.process import process
from base.utils import fix_names
import argparse
import numpy as np
from pprint import pprint
from pdb import set_trace

def collect_args():
    parser = argparse.ArgumentParser(description = "Process (prepared) JSON(s)")
    parser.add_argument('-j', '--jsons', default="all", nargs='+', type=str, help="prepared JSON(s). \"all\" will do them all. Default = all")
    parser.add_argument('-f', '--force', default=False, action='store_true', help="overwrite any intermediate files (don't restore)")
    return parser.parse_args()

def make_config (prepared_json):
    return {
        "dir": "flu",
        "in": prepared_json,
        "geo_inference": ['country', 'region'], # what traits to perform this on
        "auspice": { ## settings for auspice JSON export
            "color_options": {
                "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
                "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
            },
            "controls": {'geographic location':['country'], 'authors':['authors']}
        },
        "estimate_mutation_frequencies": [
            {"region": "global", "min_freq": 0.02, "pivot_spacing": 1.0/12, "inertia":np.exp(-1.0/12), "stiffness":0.8*12},
            {"region": "groups", "min_freq": 0.05, "inertia":np.exp(-1.0/12), "stiffness":0.8*12},
        ],
        "estimate_tree_frequencies": True,
    }

if __name__=="__main__":
    args = collect_args()
    if "all" in args.jsons:
        jsons = glob.glob("prepared/*.json")

    for prepared_json in args.jsons:
        pprint("Processing {}".format(prepared_json))
        config = make_config(prepared_json)
        if args.force:
            config["restore"] = False
        runner = process(config)
        runner.align()
        runner.estimate_mutation_frequencies_wrapper()
        runner.build_tree()
        runner.clock_filter()
        runner.annotate_tree(Tc=0.02, timetree=True, reroot='best')
        runner.run_geo_inference()

        # tree freqs could be in config like mut freqs - which is preferred?
        if config["estimate_tree_frequencies"]:
            time_interval = runner.info["time_interval"]
            pivots = np.arange(time_interval[1].year+(time_interval[1].month-1)/12.0,
                               time_interval[0].year+time_interval[0].month/12.0,
                               1.0/12.0)
            runner.estimate_tree_frequencies(pivots=pivots)
            for regionTuple in runner.info["regions"]:
                runner.estimate_tree_frequencies(region=str(regionTuple[0]))

        runner.save_as_nexus()
        runner.auspice_export()
