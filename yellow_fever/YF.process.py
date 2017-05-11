from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.process import process

config = {
    "dir": "yellow_fever",
    "output": { # will move to the default config file
        "data": "processed",
        "auspice": "auspice",
    },
    "in": "prepared/yellow-fever.json", # should be able to specify from command line
    "geo_inference": ['country'], # what traits to perform this on
    "auspice": { ## settings for auspice JSON export
        "panels": ['tree', 'map', 'entropy'],
        "extra_attr": [],
        "date_range": {'date_min': '975-01-01', 'date_max': '2017-06-01'},
        "color_options": {
            "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
            "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
            "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"}
        },
        "controls": {'geographic location':['country']}
    }
}

if __name__=="__main__":
    runner = process(config)
    runner.align()
    runner.build_tree()
    runner.clock_filter()
    runner.annotate_tree(Tc=0.02, timetree=True, reroot='best', confidence=False)
    runner.run_geo_inference()
    runner.save_as_nexus()
    runner.auspice_export()
