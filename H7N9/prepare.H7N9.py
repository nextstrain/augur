"""
This script takes fauna fasta files and produces filtered and subsetted JSONs
To be processed by augur

* RUN THIS SCRIPT IN THE H7N9 DIRECTORY!!!!

"""
from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append('reference_segments')
from reference_info import references
from base.prepare import prepare
from datetime import datetime
from base.utils import fix_names

dropped_strains = [
    "A/Chicken/Netherlands/16007311-037041/2016", # not part of epi clade
    "A/Chicken/Netherlands/1600", # not part of epi clade
    "A/duck/Zhejiang/LS02/2014", # not of part of epi clade
    "A/British Columbia/1/2015", # travel case. throws off map
    "A/BritishColumbia/1/2015", # travel case. throws off map
]

config = {
    "dir": "H7N9", # the current directory. You mush be inside this to run the script.
    "file_prefix": "flu_h7n9",
    "segments": ["pb2", "pb1", "pa", "ha", "np", "na", "mp", "ns"],
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
    # "subsample": False,
    "subsample": {
        "category": None,
        "priority": None,
        "threshold": 2,
    },
    # see the docs for what's going on with colours (sic) & lat/longs
    "colors": ["country", "division", "host"], # essential. Maybe False.
    "lat_longs": ["country", "division"], # essential. Maybe False.
    "references": references, # imported
}


if __name__=="__main__":
    runner = prepare(config)
    runner.load_references()
    runner.applyFilters()
    runner.ensure_all_segments()
    runner.subsample()
    runner.colors()
    runner.latlongs()
    runner.write_to_json()
