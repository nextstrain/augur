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

dropped_strains = [
    "ZF36_36S", # possible contamination
    "Dominican_Republic/2016/PD2", "GD01", "GDZ16001", "VEN/UF_2/2016", # true strains, but duplicates of other strains in dataset
    "Bahia04", "JAM/2016/WI_JM6", "Bahia11", "Bahia12", "DOM/2016/MA_WGS16_009", "VE_Ganxian", "BRA/2016/FC_DQ60D1", # excessive terminal branch length
    "VR10599/Pavia/2016", "34997/Pavia/2016", # exports
    "THA/PLCal_ZV/2013", "SK403/13AS", "SV0010/15", "SK364/13AS" # clock is off
]

config = {
    "dir": "zika",
    "file_prefix": "zika",
    "input_paths": ["../../fauna/data/zika.fasta"],
    "header_fields": {0:'strain', 2:'accession', 3:'date', 4:'region', 5:'country',
                    6:'division', 8:'db', 10:'authors', 11:'url'},
    "filters": (
        ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
        ("Restrict Date Range", lambda s: s.attributes['date'] >= datetime(2012,01,1).date()),
        ("Restrict Date Range", lambda s: s.attributes['date'] <= datetime(2018,01,1).date()),
        ("Sequence Length", lambda s: len(s.seq)>=2000),
    ),
    "subsample": {
        "category": lambda x:(x.attributes['date'].year, x.attributes['date'].month, x.attributes['country']),
    },
    "colors": ["country", "region"],
    "color_defs": ["./colors.tsv"],
    "lat_longs": ["country", "region"],
    "reference": {
        "path": "metadata/zika_outgroup.gb",
        "metadata": {
            'strain': "PF13/251013_18", "accession": "KX369547", "date": "2013-10-25",
            'host': "human", 'country': "French Polynesia"
        },
        "include": 2,
        "genes": ['CA', 'PRO', 'MP', 'ENV', 'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', '2K', 'NS4B', 'NS5']
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
