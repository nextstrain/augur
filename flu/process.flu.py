from __future__ import print_function
import os, sys, glob
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.process import process
from base.utils import fix_names
import argparse
from pprint import pprint


def make_config (prepared_json):
    return {
        "dir": "flu",
        "output": { # will move to the default config file
            "data": "processed",
            "auspice": "auspice",
        },
        "in": prepared_json,
        "geo_inference": ['country', 'region'], # what traits to perform this on
        "auspice": { ## settings for auspice JSON export
            "panels": ['tree', 'map', 'entropy'],
            "extra_attr": [],
            "date_range": {},
            "color_options": {
                "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
                "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
                "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
                # "ep":{"key":"ep", "legendTitle":"Epitope Mutations", "menuItem":"epitope mutations", "type":"continuous"},
                # "ne":{"key":"ne", "legendTitle":"Non-epitope Mutations", "menuItem":"nonepitope mutations", "type":"continuous"},
                # "rb":{"key":"rb", "legendTitle":"Receptor Binding Mutations", "menuItem":"RBS mutations", "type":"continuous"},
                "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"},
                # "cHI":{"key":"cHI", "legendTitle":"Antigenic advance", "menuItem":"Antigenic", "type":"continuous"}
            },
            "controls": {'geographic location':['country'], 'authors':['authors']}
        }
    }


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "Process the prepared JSONs.")
    parser.add_argument('--input_folder', default="prepared", nargs=1, type = str, help = "folder containing prepared JSONs")

    # parser.add_argument('--builds', default=['live', 'cdc'], nargs='+', type = str, help = "builds to include")
    # parser.add_argument('--lineages', default=['h3n2', 'h1n1pdm', 'vic', 'yam'], nargs='+', type = str,  help = "lineages to include")
    # parser.add_argument('--resolutions', default=['2y', '3y', '6y', '12y'], nargs='+', type = str,  help = "resolutions to include")
    # parser.add_argument('--assays', default=['hi', 'fra'], nargs='+', type = str, help = "assays to include (CDC build only)")
    # parser.add_argument('--passages', default=['egg', 'cell'], nargs='+', type = str, help = "passages to include (CDC build only)")
    # parser.add_argument('--live_auspice_path', default = '../nextflu/auspice/data/')
    # parser.add_argument('--cdc_auspice_path', default = '../nextflu-cdc/auspice/data/')
    params = parser.parse_args()

    for prepared_json in glob.glob("{}/*.json".format(params.input_folder)):
        pprint(prepared_json)
        config = make_config(prepared_json)

        runner = process(config)
        runner.align()
        runner.build_tree()
        runner.clock_filter()
        runner.annotate_tree(Tc=0.02, timetree=True, reroot='best', confidence=False)
        runner.run_geo_inference()
        runner.save_as_nexus()
        runner.auspice_export()
