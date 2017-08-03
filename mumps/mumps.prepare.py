from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.prepare import prepare
from datetime import datetime
from base.utils import fix_names
import argparse

def collect_args():
    parser = argparse.ArgumentParser(description = "Prepare fauna FASTA for analysis")
    return parser.parse_args()

dropped_strains = []

filters = {
    "dropped_strains": ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
    "canada_only": ("Canada only", lambda s: s.attributes['country'] == "canada"),
    "exclude_BC": ("Exclude BC outbreak", lambda s: not s.attributes['accession'].startswith("BC_outbreak")),
    "Mass_only": ("Massachusetts only", lambda s: s.attributes['accession'].startswith("Massachusetts_outbreak")),
    "exclude_Mass": ("Exclude Massachusetts outbreak", lambda s: not s.attributes['accession'].startswith("Massachusetts_outbreak")),
    "unknown_country": ("Exclude unknown countries", lambda s: not s.attributes['country'].startswith("unknown"))
}

def make_config(context):
    config = {
        "dir": "mumps",
        "file_prefix": "mumps_%s"%context,
        "input_paths": ["../../fauna/data/mumps.fasta"],
        "header_fields": {0:'strain', 2:'accession', 3:'date', 4:'region', 5:'country',
                        6:'division', 8:'db', 10:'authors', 11:'url'},
        "subsample": False,
        "colors": ["country", "region"],
        "color_defs": ["./colors.tsv"],
        "lat_longs": ["country", "region"],
        "reference": {
            "path": "mumps-reference.gb",
            "metadata": {
                'strain': "MuV/Gabon/13/2", "accession": "KM597072.1", "date": "2013-03-01",
                'host': "human", 'country': "Gabon"
            },
            "include": 0,
            "genes": ['NC', 'P', 'V', 'I', 'M', 'F', 'SH', 'HN', 'L']
        }
    }
    if context == "global":
        config["filters"] = (filters["dropped_strains"], filters["exclude_BC"], filters["exclude_Mass"], filters["unknown_country"])
    elif context == "bc":
        config["filters"] = (filters["dropped_strains"], filters["canada_only"], filters["unknown_country"])
    elif context == "mass":
        config["filters"] = (filters["dropped_strains"], filters["Mass_only"],filters["unknown_country"])
    else:
        print("Unknown context. FATAL")
        sys.exit(2)



    return config

if __name__=="__main__":
    params = collect_args()
    for context in ["global", "bc", "mass"]:
        runner = prepare(make_config(context))
        runner.load_references()
        runner.applyFilters()
        runner.ensure_all_segments()
        runner.subsample()
        runner.colors()
        runner.latlongs()
        runner.write_to_json()
