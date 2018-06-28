from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import base.process
from base.process import process


def collect_args():
    parser = base.process.collect_args()
    parser.set_defaults(
        json="prepared/lassa_s.json"
    )
    return parser


config = {
    "dir": "lassa",
    "subprocess_verbosity_level":3,
    "geo_inference": ['country'], # what traits to perform this on
    "auspice": { ## settings for auspice JSON export
        "panels": ['tree', 'map', 'entropy'],
        "color_options": {
            "country": {"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete", "color_map": []},
            "host_species": {"key":"host_species", "legendTitle":"Host", "menuItem":"host", "type":"discrete", "color_map": []}
        },
        "defaults": {
            "colorBy": "country",
            "geoResolution": "country"
        }
    },
    "newick_tree_options": {}, #"method":"iqtree"},
    "clock_filter":False,
    "timetree_options": {
        "Tc": "skyline",
        "resolve_polytomies": True,
        "n_points": 20,
        "stiffness": 3.0,
        "reroot":"min_dev",
        "fixed_clock_rate":0.0006
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
    runner.align()
    runner.build_tree()
    runner.timetree_setup_filter_run()
    runner.run_geo_inference()
    runner.save_as_nexus()
    runner.auspice_export()
