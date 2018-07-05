from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import base.process
from base.process import process


def collect_args():
    """Returns a WNV-specific argument parser."""
    parser = base.process.collect_args()
    return parser


config = {
    "dir": "WNV",
    "in": "prepared/WNV_NA.json",
    "newick_tree_options": {"nthreads": 3},
    "clock_filter": {
        "n_iqd": 4,
    },
    "geo_inference": ['state', "lineage"], # what traits to perform this on
    "auspice": { ## settings for auspice JSON export
        "color_options": {
            "country": {"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
            "wnv_strain": {"key":"wnv_strain", "legendTitle":"Strain", "menuItem":"strain", "type":"wnv_strain"},
            "division": {"key":"division", "legendTitle":"division", "menuItem":"division", "type":"discrete"},
            "state": {"key":"state", "legendTitle":"state", "menuItem":"state", "type":"discrete"},
            "authors": {"key":"authors", "legendTitle":"Authors", "menuItem":"authors", "type":"discrete"},
            "host": {"key":"host", "legendTitle":"Host Species", "menuItem":"host", "type":"discrete"},
        },
        "defaults": {
            'mapTriplicate': False,
            'colorBy': "state"
        }
    },
    "timetree_options": {
        "Tc": 'opt',
        "confidence":True
    }
}

if __name__=="__main__":
    parser = collect_args()
    params = parser.parse_args()

    if params.clean:
        config["clean"] = True

    if params.json:
        config["in"] = params.json

    config["newick_tree_options"]["method"] = params.tree_method

    runner = process(config)
    runner.align(fill_gaps=True)
    runner.build_tree()
    runner.timetree_setup_filter_run()
    runner.run_geo_inference()
    runner.save_as_nexus()
    runner.auspice_export()
