from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.process import process

config = {
    "dir": "mumps",
    "output": { # will move to the default config file
        "data": "processed",
        "auspice": "auspice",
    },
    "in": "prepared/mumps_SH.json", # should be able to specify from command line
    "geo_inference": ['country', 'region'], # what traits to perform this on
    "auspice": { ## settings for auspice JSON export
        "panels": ['tree', 'map', 'entropy'],
        "extra_attr": [],
        "date_range": {'date_min': '1995-01-01', 'date_max': '2017-06-01'},
        "color_options": {
            "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
            "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
            "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
            "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"},
            # "fauna_date":{"key":"fauna_date", "legendTitle":"Analysis date", "menuItem":"analysisdate", "type":"continuous"},
        },
        "controls": {'geographic location':['country'], 'authors':['authors']}
    }
}

if __name__=="__main__":
    runner = process(config)
    runner.align()
    runner.build_tree()
    runner.clock_filter()
    runner.annotate_tree()
    runner.run_geo_inference()
    runner.save_as_nexus()
    runner.auspice_export()
