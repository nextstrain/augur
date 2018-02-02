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
    parser.add_argument('-g', '--geo', default='global', type = str, help = "geo resolution, global or na")
    parser.set_defaults(
        viruses_per_month=0
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

def make_config(params):
    if params.geo == "global":
        file_prefix = "mumps_global"
        if params.viruses_per_month == 0:
            viruses_per_month = 3
        else:
            viruses_per_month = params.viruses_per_month
        dropped_strains = [
            "WA0268502_buccal/Washington.USA/16" # not yet released
        ]
        colors = ["country", "region"]
        lat_longs = ["country", "region"]
        auspice_filters = ["country", "region"]
        filters = (
            ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
            ("Restrict Date Range", lambda s: s.attributes['date'] >= datetime(1950,01,1).date()),
            ("Restrict Date Range", lambda s: s.attributes['date'] <= datetime(2020,01,1).date()),
            ("Sequence Length", lambda s: len(s.seq)>=5000),
            ("Number Ns", lambda s: s.seq.count('N')<=3000)
        )
    elif params.geo == "na":
        file_prefix = "mumps_na"
        if params.viruses_per_month == 0:
            viruses_per_month = 100
        else:
            viruses_per_month = params.viruses_per_month
        dropped_strains = [
            "Ontario.CAN/13.10/G", "Ontario.CAN/04.10/G", "Massachusetts.USA/37.16/1/G", "BritishColumbia.CAN/50.16/H",
            "BritishColumbia.CAN/22.16/1/G", "Mass.USA/4.10",
            "Virginia.USA/10.12/H", "BritishColumbia.CAN/33.16/3/G",
            "BritishColumbia.CAN/33.16/1/G", "BritishColumbia.CAN/9.17/A",
            "BritishColumbia.CAN/28.16/3/G",
            # true strains, but group outside NA outbreak clade
            "WA0268502_buccal/Washington.USA/16", # not yet released
            "Washington.USA/2017217/8.17/3/G", # outlier. MRCA with other NA strains of 1990
            "BritishColumbia.CAN/34.16/2/F" #MuV genotype F. MRCA of 1943 (!)
        ]
        colors = ["country", "division"]
        lat_longs = ["country", "division"]
        auspice_filters = ["country", "division"]
        filters = (
            ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
            ("Restrict Date Range", lambda s: s.attributes['date'] >= datetime(2009,01,1).date()),
            ("Restrict Date Range", lambda s: s.attributes['date'] <= datetime(2020,01,1).date()),
            ("Sequence Length", lambda s: len(s.seq) >= 5000),
            ("Number Ns", lambda s: s.seq.count('N') <= 3000),
            ("Restrict Region", lambda s: s.attributes['region'] == 'north_america')
        )
    config = {
        "dir": "mumps",
        "file_prefix": file_prefix,
        "title": "Real-time tracking of mumps virus evolution",
        "maintainer": ["James Hadfield", "http://bedford.io/team/james-hadfield/"],
        "input_paths": ["../../../fauna/data/mumps.fasta"],
        "header_fields": {0:'strain', 2:'accession', 3:'date', 4:'region', 5:'country',
            6:'division', 8:'db', 10:'authors', 11:'url', 12:'title',
            13: 'journal', 14: 'paper_url', 15: 'MuV_genotype'},
        "filters": filters,
        "subsample": {
            "threshold": viruses_per_month,
            "category": lambda x:(x.attributes['date'].year, x.attributes['date'].month, x.attributes['country'])
        },
        "colors": colors,
        "color_defs": ["./colors.tsv"],
        "lat_longs": lat_longs,
        "auspice_filters": auspice_filters,
        "reference": {
            "path": "mumps-reference.gb",
            "metadata": {
                'strain': "MuV/Gabon/13/2", "accession": "KM597072.1", "date": "2013-03-01",
                'host': "human", 'country': "Gabon", 'region': "Gabon", 'MuV_genotype': "G"
            },
            "include": 0,
            "genes": ['NC', 'P', 'V', 'I', 'M', 'F', 'SH', 'HN', 'L']
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
