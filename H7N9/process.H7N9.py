from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.process import process
import argparse

config = {
    "dir": "H7N9",
    "output": { # will move to the default config file
        "data": "processed",
        "auspice": "auspice",
    },
    # "in": "prepared/flu_H7N9_HA.json", # should be able to specify from command line
    "geo_inference": ['country', 'division', 'host'], # what traits to perform this on
    "geo_inference_likelihoods": True,
    "temporal_confidence": True,
    "auspice": { ## settings for auspice JSON export
        "panels": ['tree', 'map', 'entropy'],
        "extra_attr": [],
        "date_range": {},
        "analysisSlider": "fauna_date",
        "color_options": {
            "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete", "color_map": []},
            "division":{"key":"division", "legendTitle":"Division", "menuItem":"division", "type":"discrete", "color_map": []},
            "host":{"key":"host", "legendTitle":"Host", "menuItem":"host", "type":"discrete", "color_map": []},
            "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
            "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"},
            "fauna_date":{"key":"fauna_date", "legendTitle":"Analysis date", "menuItem":"Analysis date", "type":"continuous"},
        },
        "controls": {'geographic location':['country'], 'authors':['authors']}
    }
}

if __name__=="__main__":
    # for H7N9 currently things are kinda easy - everything needed is defined in the JSON!
    parser = argparse.ArgumentParser(description = "Process a given JSONs")
    parser.add_argument('-j', '--json', required=True, type=str, help="prepared JSON")
    params = parser.parse_args()

    config["in"] = params.json
    runner = process(config)
    runner.align()
    runner.build_tree()
    runner.clock_filter()
    runner.annotate_tree(Tc=0.02, timetree=True, reroot='best')
    runner.run_geo_inference()
    runner.save_as_nexus()
    runner.auspice_export()
