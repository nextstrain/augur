from __future__ import division, print_function
from pprint import pprint

"""
Here are the default config dictionaries.
The config dictionaries provided in each pathogen's prepare / process
file (e.g. process.H7N9.py) overwrite this.
This allows for smaller / easier to interpret scripts.
"""

prepare = {
    "segments": False,
    "input_format": "fasta",
    "output_folder": "prepared",
    "date_format": ["%Y-%m-%d"],
    "require_dates": True,
    "subsample": False,
}

process = {
    "output": {
        "data": "processed",
        "auspice": "auspice",
    },
    "geo_inference_likelihoods": True,
    "temporal_confidence": True,
    "auspice": { ## settings for auspice JSON export
        "panels": ['tree', 'map', 'entropy'],
        "extra_attr": [],
        "color_options": {
            "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
            "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"}
        },
    }
}

def merge(a, b):
    for key, value in b.iteritems():
        if key in a and isinstance(value, dict) and isinstance(b[key], dict):
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
    # they can be merged together...
    config = merge(config, user_config)
    # pprint(config)
    # pprint(config["auspice"])
    return config
