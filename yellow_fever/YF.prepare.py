from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.prepare import prepare
from datetime import datetime
from base.utils import fix_names

dropped_strains = []

config = {
    "dir": "yellow_fever",
    "file_prefix": "yellow-fever",
    "input_paths": ["../../fauna/data/yellow_fever.fasta"],
    # ArD114972|JX898872|1995-XX-XX|senegal|mosquito
    "header_fields": {0:'strain', 1:'accession', 2:'date', 3:'country', 4:'host'},
    "filters": (
        ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
    ),
    "colors": ["country"], # essential. Maybe False.
    "lat_longs": ["country"], # essential. Maybe False.
    "lat_long_defs": '../../fauna/source-data/geo_lat_long.tsv',
    "reference": {
        "path": "YF.reference.gb",
        "metadata": {
            'strain': "reference", "accession": "NC_002031", "date": "XXXX-XX-XX",
            'host': "?", 'country': "?"
        },
        "use": False,
        "genes": ['poly']
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
