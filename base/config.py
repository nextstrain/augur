from __future__ import division, print_function
from pprint import pprint

"""
Here are the default config dictionaries.
The config dictionaries provided in each pathogen's prepare / process
file (e.g. process.H7N9.py) overwrite this.
This allows for smaller / easier to interpret scripts.
"""

required_config_fields = {
    "prepare": [
        "dir", "file_prefix", "segments", "input_format", "input_paths",
        "output_folder", "header_fields", "date_format",
        "require_dates", "colors", "lat_longs"
    ],
    "process": [
    ]
}


prepare = {
    "segments": False,
    "input_format": "fasta",
    "output_folder": "prepared",
    "date_format": ["%Y-%m-%d"],
    "require_dates": True,
    "subsample": False,
    "ensure_all_segments": True, #this is ignored if only 1 segment
    "lat_long_defs": '../../../fauna/source-data/geo_lat_long.tsv',
    "maintainer": "unknown"
}

process = {
    "output": {
        "data": "processed",
        "auspice": "auspice",
    },
    "geo_inference": False,
    "geo_inference_options": {
        "confidence": True
    },
    "temporal_confidence": True,
    "auspice": { ## settings for auspice JSON export
        "panels": ['tree', 'map', 'entropy'],
        "extra_attr": [],
        "color_options": {
            "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
            "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"}
        },
        "controls": {}
    },
    "clock_filter": {
        "n_iqd": 3,
        "plot": False,
        "remove_deep_splits": False,
    },
    "estimate_mutation_frequencies": False,
    "estimate_tree_frequencies": False,
    "titers": False,
    "subprocess_verbosity_level": 0,
    "newick_tree_options": {},
    "timetree_options": {
        "Tc": 0.02,
        "reroot": 'best',
    },
    "clean": False,
}

def merge(a, b):
    for key, value in b.iteritems():
        if key in a and isinstance(a[key], dict) and isinstance(b[key], dict):
            a[key] = merge(a[key], b[key])
        else:
            a[key] = b[key]
    return a


def combine_configs(config_type, user_config):
    if config_type == "prepare":
        config = prepare.copy()
    elif config_type == "process":
        config = process.copy()
    else:
        raise Exception("unknown config type demanded: {}".format(config_type))

    # instead of just doing
    # config.update(user_config)
    # update key by key so that if the value is itself a dictionary
    # they can be merged together... (recursively)
    config = merge(config, user_config)
    if "geo_inference" in config and config["geo_inference"] != False and "geographic location" not in config["auspice"]["controls"]:
        config["auspice"]["controls"]["geographic location"] = config["geo_inference"]


    if config_type == "prepare" and "title" not in config:
        config["title"] = config["file_prefix"]

    # pprint(config)
    # pprint(config["auspice"])

    # check config file appears to be OK
    try:
        for x in required_config_fields[config_type]:
            assert(x in config)
    except AssertionError:
        print("Fatal Error: Config file is missing field '{}'".format(x))
        sys.exit(2)
    return config
