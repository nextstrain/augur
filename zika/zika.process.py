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
    "in": "prepared/zika.json", # should be able to specify from command line
    "geo_inference": ['country', 'region'], # what traits to perform this on
    "auspice": { ## settings for auspice JSON export
        "color_options": {
            "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
            "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
        },
        "controls": {'authors':['authors']}
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
    runner.annotate_tree()
    runner.run_geo_inference()
    runner.save_as_nexus()
    runner.auspice_export()
