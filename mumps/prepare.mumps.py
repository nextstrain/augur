"""
This script takes fauna fasta files and produces filtered and subsetted JSONs
To be processed by augur

* RUN THIS SCRIPT IN THE MUMPS DIRECTORY!!!!
"""
from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.prepare import prepare
from datetime import datetime
from base.utils import fix_names

dropped_strains = ["Jeryl_Lynn_Vaccine_Strain_X63707"]

config = {
    "dir": "mumps", # the current directory. You mush be inside this to run the script.
    "file_prefix": "mumps",
    "segments": ["SH"],
    "input_paths": ["test_input/284_SH_Jenn.mfa"],
    "header_fields": {0:'strain', 1:'region', 2:'country', 3:'date'},
    "date_format": ["%Y-%W"],
    "filters": (
        ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
    ),
    "subsample": False,
    "colors": ["country", "region"],
    "lat_longs": ["country", "region"],
    "lat_long_defs": ['../../fauna/source-data/geo_lat_long.tsv', "geo_lat_long.tsv"],
    "references": {
        "SH": {
            "path": "reference.gb",
            "metadata": {
                'strain': "MuVs_KF840226", "date": "2012-21",
                'region': "Goteborg", 'country': "SWE"
            },
            "include": 2,
            "genes": ["SH"]
        }
    }
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
