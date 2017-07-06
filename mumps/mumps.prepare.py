from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.prepare import prepare
from datetime import datetime
from base.utils import fix_names
import argparse

def collect_args():
    parser = argparse.ArgumentParser(description = "Prepare fauna FASTA for analysis")
    parser.add_argument('-v', '--viruses_per_month', type = int, default = 15, help='Subsample x viruses per country per month. Set to 0 to disable subsampling. (default: 15)')
    return parser.parse_args()

dropped_strains = []

config = {
    "dir": "mumps",
    "file_prefix": "mumps",
    "input_paths": ["../../fauna/data/mumps.fasta"],
    "header_fields": {0:'strain', 2:'accession', 3:'date', 4:'region', 5:'country',
                    6:'division', 8:'db', 10:'authors', 11:'url'},
    "filters": (
        ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
    ),
    "subsample": {
        "category": lambda x:(x.attributes['date'].year, x.attributes['date'].month),
    },
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

if __name__=="__main__":
    params = collect_args()
    if params.viruses_per_month == 0:
        config["subsample"] = False
    else:
        config["subsample"]["threshold"] = params.viruses_per_month
    runner = prepare(config)
    runner.load_references()
    runner.applyFilters()
    runner.ensure_all_segments()
    runner.subsample()
    runner.colors()
    runner.latlongs()
    runner.write_to_json()
