from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.process import process
import argparse

parser = argparse.ArgumentParser(description = "Process a given JSONs")
parser.add_argument('-r', '--restore', action='store_true', help="try to restore")
params = parser.parse_args()


config = {
    "dir": "zika",
    "output": { # will move to the default config file
        "data": "processed",
        "auspice": "auspice",
    },
    "in": "prepared/zika.json", # should be able to specify from command line
    "geo_inference": ['country', 'region'], # what traits to perform this on
    "geo_inference_likelihoods": True,
    "temporal_confidence": True,
    "auspice": { ## settings for auspice JSON export
        "panels": ['tree', 'map', 'entropy'],
        "extra_attr": [],
        "date_range": {'date_min': '2013-06-01', 'date_max': '2017-06-01'},
        "color_options": {
            "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
            "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
            "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
            "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"}
        },
        "controls": {'geographic location':['region', 'country'], 'authors':['authors']}
    }
}

if __name__=="__main__":
    params = parser.parse_args()
    runner = process(config)
    runner.align()
    # if not params.restore:
    #     runner.align()
    #     runner.dump()
    # else:
    #     runner.load()
    runner.build_tree()
    runner.clock_filter()
    runner.annotate_tree(Tc=0.02, timetree=True, reroot='best')
    runner.run_geo_inference()
    runner.save_as_nexus()
    runner.auspice_export()
