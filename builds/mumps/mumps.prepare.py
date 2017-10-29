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
    parser.set_defaults(
        viruses_per_month=15,
        file_prefix="mumps"
    )
    return parser.parse_args()

# filters = {
#     "dropped_strains": ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
#     "canada_only": ("Canada only", lambda s: s.attributes['country'] == "canada"),
#     "exclude_BC": ("Exclude BC outbreak", lambda s: not s.attributes['accession'].startswith("BC_outbreak")),
#     "Mass_only": ("Massachusetts only", lambda s: s.attributes['accession'].startswith("Massachusetts_outbreak")),
#     "exclude_Mass": ("Exclude Massachusetts outbreak", lambda s: not s.attributes['accession'].startswith("Massachusetts_outbreak")),
#     "unknown_country": ("Exclude unknown countries", lambda s: not s.attributes['country'].startswith("unknown"))
# }

dropped_strains = [
]

def make_config():
    config = {
        "dir": "mumps",
        "file_prefix": "mumps",
        "title": "Real-time tracking of mumps virus evolution",
        "maintainer": ["James Hadfield", "http://bedford.io/team/james-hadfield/"],
        "input_paths": ["../../../fauna/data/mumps.fasta"],
        "header_fields": {0:'strain', 2:'accession', 3:'date', 4:'region', 5:'country',
            6:'division', 8:'db', 10:'authors', 11:'url', 12:'title',
            13: 'journal', 14: 'paper_url'},
        "filters": (
            ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
            ("Restrict Date Range", lambda s: s.attributes['date'] >= datetime(1950,01,1).date()),
            ("Restrict Date Range", lambda s: s.attributes['date'] <= datetime(2020,01,1).date()),
            ("Sequence Length", lambda s: len(s.seq)>=5000),
            ("Number Ns", lambda s: s.seq.count('N')<=3000)
        ),
        "subsample": {
            "category": lambda x:(x.attributes['date'].year, x.attributes['date'].month, x.attributes['country']),
        },
        "colors": ["country", "region"],
        "color_defs": ["./colors.tsv"],
        "lat_longs": ["country", "region"],
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
