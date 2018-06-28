from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import base.process
from base.process import process
import argparse


def collect_args():
    """Returns an avian flu-specific argument parser.
    """
    parser = base.process.collect_args()
    parser.set_defaults(
        json="prepared/flu_avian_h7n9_ha.json"
    )
    return parser


config = {
    "dir": "avian",
    # "in": "prepared/flu_H7N9_HA.json", # should be able to specify from command line
    "geo_inference": ["division"], # what traits to perform this on
    "auspice": { ## settings for auspice JSON export
        "panels": ['tree', 'map', 'entropy'],
        "defaults": {
            "colorBy": "division",
            "geoResolution": "division"
        },
        "extra_attr": [],
        "date_range": {},
        "color_options": {
            "division":{"key":"division", "legendTitle":"Division", "menuItem":"division", "type":"discrete", "color_map": []},
            "host":{"key":"host", "legendTitle":"Host", "menuItem":"host", "type":"discrete", "color_map": []},
            "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
            "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"}
        }
    },
    "newick_tree_options": {}
}

if __name__=="__main__":
    parser = collect_args()
    params = parser.parse_args()

    config["in"] = params.json
    config["clean"] = params.clean

    runner = process(config)
    runner.align()
    runner.build_tree()
    runner.timetree_setup_filter_run()
    runner.run_geo_inference()
    runner.save_as_nexus()
    runner.auspice_export()
