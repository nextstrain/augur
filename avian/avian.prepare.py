"""
This script takes fauna fasta files and produces filtered and subsetted JSONs
To be processed by augur

* RUN THIS SCRIPT IN THE avian DIRECTORY!!!!

"""
from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside avian or zika
sys.path.append('reference_segments')
from reference_info import references
import base.prepare
from base.prepare import prepare
from datetime import datetime
from base.utils import fix_names
import argparse


def collect_args():
    """Returns an avian flu-specific argument parser.
    """
    parser = base.prepare.collect_args()

    segments = ["pb2", "pb1", "pa", "ha", "np", "na", "mp", "ns"]
    parser.add_argument('-l', '--lineage', choices=['h7n9'], default='h7n9', type=str, help="serotype")
    parser.add_argument('-s', '--segments', choices=segments, default=segments, nargs='+', type = str,  help = "segments to prepare")
    parser.set_defaults(
        viruses_per_month=0,
        file_prefix="avian_h7n9"
    )

    return parser

dropped_strains = [
    "A/Chicken/Netherlands/16007311-037041/2016", # not part of epi clade
    "A/Chicken/Netherlands/1600", # not part of epi clade
    "A/duck/Zhejiang/LS02/2014", # not of part of epi clade
    "A/British Columbia/1/2015", # travel case. throws off map
    "A/BritishColumbia/1/2015", # travel case. throws off map
]

config = {
    "dir": "avian", # the current directory. You must be inside this to run the script.
    "file_prefix": "avian_h7n9",
    "title": "Genomic analysis of the 5 Epidemics Influenza A/H7N9",
    "maintainer": ["@jh_viz", "https://twitter.com/jh_viz"],
    "input_paths": [
        "../../fauna/data/h7n9_pb2.fasta",
        "../../fauna/data/h7n9_pb1.fasta",
        "../../fauna/data/h7n9_pa.fasta",
        "../../fauna/data/h7n9_ha.fasta",
        "../../fauna/data/h7n9_np.fasta",
        "../../fauna/data/h7n9_na.fasta",
        "../../fauna/data/h7n9_mp.fasta",
        "../../fauna/data/h7n9_ns.fasta"
    ],
    "header_fields": {
        0:'strain', 2:'accession', 3:'date', 4:'host', 6:'country', 7: 'division', 11: 'fauna_date'
    },
    "traits_are_dates": ["fauna_date"],
    "filters": (
        ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
        ("Prior to Epidemic", lambda s: s.attributes['date'] >= datetime(2013,1,1).date()),
        ("Missing Month Data", lambda s: "-XX-XX" not in s.attributes['raw_date']),
        ("Exclude bad host", lambda s: s.attributes["host"] not in ["laboratoryderived", "watersample"]),
        # ("Restrict to Humans", lambda s: s.attributes["host"] in ["human"]),
        ("Sequence Length", {
            "pb2": lambda s: len(s.seq)>=2200,
            "pb1": lambda s: len(s.seq)>=2200,
            "pa": lambda s: len(s.seq)>=2100,
            "ha": lambda s: len(s.seq)>=1500,
            "np": lambda s: len(s.seq)>=1400,
            "na": lambda s: len(s.seq)>=1200,
            "mp": lambda s: len(s.seq)>=900,
            "ns": lambda s: len(s.seq)>=800
        })
    ),
    "subsample": {
        "category": lambda x:(x.attributes['date'].year, x.attributes['date'].month, x.attributes['division']),
        "priority": None,
    },
    # see the docs for what's going on with colours (sic) & lat/longs
    "colors": ["country", "division", "host"], # essential. Maybe False.
    "lat_longs": ["country", "division"], # essential. Maybe False.
    "references": references, # imported
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

    config["segments"] = params.segments

    runner = prepare(config)
    runner.load_references()
    runner.applyFilters()
    runner.ensure_all_segments()
    runner.subsample()
    runner.colors()
    runner.latlongs()
    runner.write_to_json()
