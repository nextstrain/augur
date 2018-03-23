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
    """Returns a Measles-specific argument parser.
    """
    parser = base.prepare.collect_args()
    parser.set_defaults(
        viruses_per_month=10
    )
    return parser.parse_args()

def make_config(params):
    dropped_strains = [
        "temara.MOR/24.03", "Mvs/Toulon.FRA/08.07" # clock is off
    ]
    filters = (
        ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
        ("Restrict Date Range", lambda s: s.attributes['date'] >= datetime(1950,01,1).date()),
        ("Restrict Date Range", lambda s: s.attributes['date'] <= datetime(2020,01,1).date()),
        ("Sequence Length", lambda s: len(s.seq)>=5000),
        ("Number Ns", lambda s: s.seq.count('N')<=3000)
    )
    config = {
        "dir": "measles",
        "file_prefix": "measles",
        "title": "Real-time tracking of measles virus evolution",
        "maintainer": ["Trevor Bedford", "http://bedford.io/team/trevor-bedford/"],
        "input_paths": ["../../../fauna/data/measles.fasta"],
        "header_fields": {0:'strain', 2:'accession', 3:'date', 4:'region', 5:'country',
            6:'division', 8:'db', 10:'authors', 11:'url', 12:'title',
            13: 'journal', 14: 'paper_url'},
        "filters": filters,
        "subsample": {
            "threshold": params.viruses_per_month,
            "category": lambda x:(x.attributes['date'].year, x.attributes['date'].month, x.attributes['country'])
        },
        "colors": ["authors", "country", "region"],
        "color_defs": ["./colors.tsv"],
        "lat_longs": ["country", "region"],
        "auspice_filters": ["authors", "region", "country"],
        "reference": {
            "path": "measles-reference.gb",
            "metadata": {
                'strain': "Ichinose-B95a", "accession": "NC_001498.1", "date": "XXXX-XX-XX",
                'host': "human", 'country': "Unknown", 'region': "Unknown"
            },
            "include": 0,
            "genes": ['N', 'P', 'V', 'C', 'M', 'F', 'H', 'L']
        }
    }

    return config

if __name__=="__main__":
    params = collect_args()
    runner = prepare(make_config(params))
    runner.load_references()
    runner.applyFilters()
    runner.ensure_all_segments()
    runner.subsample()
    runner.colors()
    runner.latlongs()
    runner.write_to_json()
