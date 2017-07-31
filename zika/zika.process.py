from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
import base.process
from base.process import process
import argparse


def collect_args():
    parser = base.process.collect_args()
    parser.set_defaults(
        json="prepared/zika.json"
    )
    return parser


config = {
    "dir": "zika",
    "in": "prepared/zika.json",
    "newick_tree_options": {"nthreads": 4},
    "clock_filter": {
        "n_iqd": 4,
    },
    "geo_inference": ['country', 'region'], # what traits to perform this on
    "geo_inference_options": {
        "root_state": {
            "region": "southeast_asia",
            "country": "vietnam",
        },
    },
    "auspice": { ## settings for auspice JSON export
        "color_options": {
            "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
            "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
        },
        "controls": {'authors':['authors']},
        "defaults": {'mapTriplicate': True},
    },
    "timetree_options": {
        "Tc": 'opt',
        "confidence":True,
    },
    "newick_tree_options":{
        "raxml": True
    }
}

if __name__=="__main__":
    parser = collect_args()
    params = parser.parse_args()

    if params.clean:
        config["clean"] = True

    if params.json:
        config["in"] = params.json

    runner = process(config)
    runner.align(fill_gaps=True)
    runner.build_tree()
    runner.timetree_setup_filter_run()
    runner.run_geo_inference()
    runner.save_as_nexus()
    runner.auspice_export()
