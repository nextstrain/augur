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
    """Returns a Zika-specific argument parser.
    """
    parser = base.prepare.collect_args()
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

def make_config():
    config = {
        "dir": "mumps",
        "file_prefix": "mumps",
        "title": "Mumps virus",
        "maintainer": ["@jh_viz", "https://twitter.com/jh_viz"],
        "input_paths": ["../../../fauna/data/mumps.fasta"],

        "header_fields": {
            0: 'strain',
            2: 'accession',
            3: 'date',
            4: 'country',
            5: 'region',
            6: 'muv_genotype',
            7: 'host',
            8: 'authors',
            9: 'title',
            10: 'journal',
            11: 'puburl',
            12: 'url'
        },
        "subsample": False,
        "colors": ["country", "region"],
        "lat_longs": ["country", "region"],
        "lat_long_defs": './geo_lat_long.tsv',
        "filters": (
            ("Sequence Length", lambda s: len(s.seq)>=13000),
            ("number Ns", lambda s: s.seq.count('N')<=3000)
        ),
        "auspice_filters": ["country", "region"],
        "reference": {
            "path": "mumps-reference.gb",
            "metadata": {
                'strain': "MuV/Gabon/13/2", "accession": "KM597072.1", "date": "2013-03-01",
                'host': "human", 'country': "Gabon", 'region': "Gabon"
            },
            "include": 0,
            "genes": ['NC', 'P', 'V', 'I', 'M', 'F', 'SH', 'HN', 'L']
        }
    }


    return config

if __name__=="__main__":
    params = collect_args()
    runner = prepare(make_config())
    runner.load_references()
    runner.applyFilters()
    runner.ensure_all_segments()
    runner.subsample()
    runner.colors()
    runner.latlongs()
    runner.write_to_json()
