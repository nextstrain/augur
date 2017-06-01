from __future__ import print_function
import os, sys, glob
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.process import process
from base.utils import fix_names
import argparse
from pprint import pprint
from pdb import set_trace

def collect_args():
    parser = argparse.ArgumentParser(description = "Process (prepared) JSON(s)")
    parser.add_argument('-j', '--jsons', default="all", nargs='+', type=str, help="prepared JSON(s). \"all\" will do them all. Default = all")
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
        }
    }

if __name__=="__main__":
    jsons = collect_args().jsons
    if "all" in jsons:
        jsons = glob.glob("prepared/*.json")

    for prepared_json in jsons:
        pprint("Processing {}".format(prepared_json))
        config = make_config(prepared_json)

        runner = process(config)
        runner.align()
        runner.build_tree()
        runner.clock_filter()
        runner.annotate_tree(Tc=0.02, timetree=True, reroot='best')
        runner.run_geo_inference()
        runner.save_as_nexus()
        runner.auspice_export()
