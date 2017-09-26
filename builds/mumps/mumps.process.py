from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
from base.process import process
import argparse
import glob

parser = argparse.ArgumentParser(description = "Process a prepared zika JSON")
parser.add_argument('--clean', action='store_true', help="clean build (remove previous checkpoints)")
parser.add_argument('-j', '--jsons', default=["all"], nargs='+', type=str, help="prepared JSON(s). \"all\" will do them all. Default = all")

def make_config(prepared_json, clean):
    return {
        "dir": "mumps",
        "in": prepared_json,
        "newick_tree_options": {"nthreads": 4},
        "clock_filter": {
            "n_iqd": 4,
        },
        "geo_inference": False,
        "auspice": { ## settings for auspice JSON export
            "color_options": {
                "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
                "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
            },
        },
        "clean": clean,
    }

if __name__=="__main__":
    params = parser.parse_args()
    jsons = glob.glob("prepared/*.json") if "all" in params.jsons else params.jsons
    for prepared_json in jsons:
        print("Processing {}".format(prepared_json))
        runner = process(make_config(prepared_json, params.clean))
        runner.align()
        runner.build_tree()
        runner.timetree_setup_filter_run()
        runner.run_geo_inference()
        runner.save_as_nexus()
        runner.auspice_export()
