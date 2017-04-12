"""
This script takes fauna fasta files and produces filtered and subsetted JSONs
To be processed by augur

* RUN THIS SCRIPT IN THE H7N9 DIRECTORY!!!!
* The reference sequence will _never_ be filtered / subsampled e.t.c.
* However, during processing, you can drop the reference for analysis

"""
from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.prepare import prepare
from datetime import datetime

dropped_strains = []

config = {
    "dir": "H7N9", # the current directory. You mush be inside this to run the script.
    "file_prefix": "H7N9",
    "segments": ["HA", "NA"],
    # "segments": ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"], # set to False, or ["something"] if not segmented...
    "input_format": "fasta",
    "input_paths": [
        # "../../fauna/data/h7n9_pb2.fasta",
        # "../../fauna/data/h7n9_pb1.fasta",
        # "../../fauna/data/h7n9_pa.fasta",
        # "../../fauna/data/h7n9_ha.fasta",
        # "../../fauna/data/h7n9_np.fasta",
        # "../../fauna/data/h7n9_na.fasta",
        # "../../fauna/data/h7n9_mp.fasta",
        # "../../fauna/data/h7n9_ns.fasta",
        "test_input/h7n9.ha.fasta",
        "test_input/h7n9.na.fasta"
    ],
    "output_folder": "prepared",
    # note that "strain" is essential and "date" is special
    "header_fields": {
        0:'strain', 2:'accession', 3:'date', 4:'host', 6:'country', 7: 'division', 11: 'fauna_date'
    },
    # can provide multiple date formats here - FIFO
    "date_format": ["%Y-%m-%d"],
    "require_dates": True,

    # require that all selected isolates have all the genomes
    "ensure_all_segments": True,

    # see docs for help with filters - nested tuples abound
    "filters": (
        ("Dropped Strains", lambda s: s.id not in [
            "A/Chicken/Netherlands/16007311-037041/2016", # not part of epi clade
            "A/Chicken/Netherlands/1600", # not part of epi clade
            "A/duck/Zhejiang/LS02/2014", # not of part of epi clade
            "A/BritishColumbia/1/2015" # travel case. throws off map
        ]),
        ("Prior to Epidemic", lambda s: s.attributes['date'] >= datetime(2013,1,1).date()),
        ("Exclude bad host", lambda s: s.attributes["host"] not in ["laboratoryderived", "watersample"]),
        ("Sequence Length", {
            "HA": lambda s: len(s.seq)>=1500,
            "NA": lambda s: len(s.seq)>=1200
        })
    ),
    # see the docs for this too! if you don't want to subsample, set it to False
    # "subsample": False,
    "subsample": {
        "category": None,
        "priority": None,
        "threshold": None,
    },

    # see the docs for what's going on with colours (sic) & lat/longs
    "colors": ["country", "division"]
}


if __name__=="__main__":
    runner = prepare(config)
    runner.applyFilters()
    runner.ensure_all_segments()
    runner.subsample()
    runner.colors()
    runner.latlongs()
    runner.write_to_json()
