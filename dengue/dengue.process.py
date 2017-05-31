from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.process import process
import argparse

def collect_args():
    parser = argparse.ArgumentParser(description = "Prepare fauna FASTA for analysis")
    parser.add_argument('-j', '--json', required=True, type=str, help="prepared JSON file")
    return parser.parse_args()

config = {
    "dir": "dengue",
    "geo_inference": ['region'],
    "auspice": {
        "color_options": {
            "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
        },
        "controls": {'geographic location':['region'], 'authors':['authors']},
        "defaults": {'mapTriplicate': True, 'geoResolution': 'region'}
    }
}

if __name__=="__main__":
    params = collect_args()
    config["in"] = params.json
    runner = process(config)
    runner.align()
    runner.build_tree()
    runner.clock_filter()
    runner.annotate_tree(Tc=0.02, timetree=True, reroot='best')
    runner.run_geo_inference()
    runner.save_as_nexus()
    runner.auspice_export()
