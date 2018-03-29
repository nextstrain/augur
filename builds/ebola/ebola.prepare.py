from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import base.prepare
from base.prepare import prepare
from datetime import datetime
from base.utils import fix_names
import argparse

dropped_strains = [
    "EM_074335", "EM_074462" # lack country metadata
]
forced_strains = [
    "EM_076610" # flare-up index case
]

def collect_args():
    """Returns an Ebola-specific argument parser.
    """
    parser = base.prepare.collect_args()
    parser.set_defaults(
        viruses_per_month=100,
        file_prefix="ebola"
    )
    return parser

def make_config(params):
    return {
        "dir": "ebola",
        "file_prefix": "ebola",
        "title": "Genomic epidemiology of the 2013-2016 West African Ebola epidemic",
        "maintainer": ["Trevor Bedford", "http://bedford.io/team/trevor-bedford/"],
        "input_paths": ["../../../fauna/data/ebola.fasta"],
        "header_fields": {0:'strain', 2:'accession', 3:'date', 4:'region', 5:'country', 6:'division', 8:'db', 10:'authors', 11:'url'},
        "filters": (
            ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
            ("Restrict Date Range", lambda s: s.attributes['date'] >= datetime(2012,01,1).date()),
            ("Restrict Date Range", lambda s: s.attributes['date'] <= datetime(2018,01,1).date()),
        ),
        "subsample": {
            "category": lambda x:(x.attributes['region'], x.attributes['date'].year, x.attributes['date'].month),
            "threshold": params.viruses_per_month,
            "priority": lambda x:x.id in forced_strains
        },
        "colors": ["authors", "country", "division"], # essential. Maybe False.
        "color_defs": ["./colors.tsv"],
        "lat_longs": ["country", "division"], # essential. Maybe False.
        "auspice_filters": ["authors", "country", "division"],
        "reference": {
            "path": "metadata/ebola_outgroup.gb",
            "metadata": {
                'strain': "reference", "accession": "KR075003", "date": "2014-XX-XX",
                'host': "human", 'country': "Liberia"
            },
            "include": 0,
            "genes": ['NP', 'VP35', 'VP40', 'GP', 'sGP', 'VP30', 'VP24', 'L']
        }
    }


if __name__=="__main__":
    parser = collect_args()
    params = parser.parse_args()
    config = make_config(params)

    if params.viruses_per_month == 0:
        config["subsample"] = False

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
