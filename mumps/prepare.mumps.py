"""
This script takes fauna fasta files and produces filtered and subsetted JSONs
To be processed by augur

* RUN THIS SCRIPT IN THE H7N9 DIRECTORY!!!!

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
    "segments": ["SH"], # set to False, or ["something"] if not segmented...
    "input_format": "fasta",
    "input_paths": ["test_input/284_SH_Jenn.mfa"],
    "output_folder": "prepared",
    "header_fields": {0:'strain', 1:'region', 2:'country', 3:'date'},
    "date_format": ["%Y-%W"],
    "require_dates": True,

    # require that all selected isolates have all the genomes
    "ensure_all_segments": True,

    # see docs for help with filters - nested tuples abound
    "filters": (
        ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
    ),
    # see the docs for this too! if you don't want to subsample, set it to False
    "subsample": False,
    # "subsample": {
    #     "category": None,
    #     "priority": None,
    #     "threshold": 1,
    # },

    # see the docs for what's going on with colours (sic) & lat/longs
    "colors": ["country", "region"], # essential. Maybe False.
    "lat_longs": ["country", "region"], # essential. Maybe False.
    "lat_long_defs": ['../../fauna/source-data/geo_lat_long.tsv',"geo_lat_long.tsv"],

    # again, see docs for reference definitions
    "references": {
        "SH": {
            "path": "reference.gb",
            "metadata": {
                'strain': "MuVs_KF840226", "date": "2012-21",
                'region': "Goteborg", 'country': "SWE"
            },
            "use": True,
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
