from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.process import process
import argparse

config = {
    "dir": "ebola",
    "in": "prepared/ebola.json",
    "geo_inference": ['country', 'division'], # what traits to perform this on
    "auspice": { ## settings for auspice JSON export
        "color_options": {
            "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
            "division":{"key":"division", "legendTitle":"Division", "menuItem":"division", "type":"discrete"},
        },
        "controls": {'geographic location':['country', 'division'], 'authors':['authors']},
        "defaults": {'mapTriplicate': False}
    },
    "timetree_options": {
        "Tc": "skyline",
        "resolve_polytomies": True,
        "n_points": 20,
        "stiffness": 3.0,
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
