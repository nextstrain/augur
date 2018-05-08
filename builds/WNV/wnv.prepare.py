from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import base.prepare
from base.prepare import prepare
from datetime import datetime
from base.utils import fix_names
import argparse

def collect_args():
    """Returns a WNV-specific argument parser."""
    parser = base.prepare.collect_args()
    parser.set_defaults(
        viruses_per_month=1,
    )
    return parser

dropped_strains = [
]

config = {
    "dir": "WNV",
    "file_prefix": "WNV_NA",
    "title": "Twenty years of West Nile Virus in North America",
    "maintainer": ["James Hadfield / Nate Grubaugh", ""],
    "input_paths": ["./data/WNV.fasta"],
    # >W112|2016-02-18|USA|CA|SD|-117.0239963|32.6299784
    "header_fields": {0:'strain', 1:'date', 2:'country', 3:'state', 4:'division', 7:'host', 8: 'wnv_strain', 9: 'authors', 10:'journal', 11:'title', 12:'url'},
    "filters": (
        ("Remove Shabman et al seqs (JCVI data!)", lambda s: s.attributes["authors"] != "Shabman_et_al"),
        ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
        # ("Restrict Date Range", lambda s: s.attributes['date'] >= datetime(2012,01,1).date()),
        # ("Restrict Date Range", lambda s: s.attributes['date'] <= datetime(2018,01,1).date()),
        ("Sequence Length", lambda s: len(s.seq)>=10000),
    ),
    "subsample": {
        "category": lambda x:(x.attributes['date'].year, x.attributes['date'].month),
    },
    "colors": ["country", "division", "state", "authors", "host", "wnv_strain"],
    "color_defs": ["./colors.tsv"],
    "lat_longs": ["state"],
    "lat_long_defs": "./lat_longs.tsv",
    "auspice_filters": ["country", "division", "state", "authors", "wnv_strain"],
    "reference": {
        "path": "reference.gb",
        "metadata": {
            'strain': "NY99", "date": "1999-XX-XX", "country": "USA", "state": "NY", "division": "?"
        },
        "include": 0,
        "genes": ['capsid', 'prM', 'env', 'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', '2K', 'NS4B', "NS5"]
    }
}

if __name__=="__main__":
    parser = collect_args()
    params = parser.parse_args()
    if params.viruses_per_month == 0:
        config["subsample"] = False
    else:
        config["subsample"]["threshold"] = params.viruses_per_month

    if params.sequences is not None:
        config["input_paths"] = params.sequences

    if params.file_prefix is not None:
        config["file_prefix"] = params.file_prefix

    runner = prepare(config)
    runner.load_references()
    runner.applyFilters()
    runner.ensure_all_segments()
    runner.subsample()
    runner.colors()
    runner.latlongs()
    runner.write_to_json()
