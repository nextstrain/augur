from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.prepare import prepare
from datetime import datetime
from base.utils import fix_names

dropped_strains = [
    "THA/PLCal_ZV/2013", "PLCal_ZV", # true strains, too basal for analysis
    "ZF36_36S", # possible contamination
    "Dominican_Republic/2016/PD2", "GD01", "GDZ16001", "VEN/UF_2/2016", # true strains, but duplicates of other strains in dataset
    "Bahia04", "JAM/2016/WI_JM6", # excessive terminal branch length
    "THA/2014/SV0127_14", "ZK_YN001", "NIID123/2016", # true strains, too basal for analysis
    "ZKA_16_291", "ZKA_16_097" # singapore, too basal for analysis
]

config = {
    "dir": "zika",
    "file_prefix": "zika",
    "segments": False,
    "input_format": "fasta",
    "input_paths": ["../../fauna/data/zika.fasta"],
    "output_folder": "prepared",
    "header_fields": {0:'strain', 2:'accession', 3:'date', 4:'region', 5:'country',
                    6:'division', 8:'db', 10:'authors', 11:'url'},
    "date_format": ["%Y-%m-%d"],
    "require_dates": True,
    "filters": (
        ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
        ("Restrict Date Range", lambda s: s.attributes['date'] >= datetime(2012,01,1).date()),
        ("Restrict Date Range", lambda s: s.attributes['date'] <= datetime(2018,01,1).date()),
        ("Sequence Length", lambda s: len(s.seq)>=2000),
    ),
    "subsample": {
        "category": lambda x:(x.attributes['date'].year, x.attributes['date'].month),
        "threshold": 100,
    },
    "colors": ["country", "division"], # essential. Maybe False.
    "lat_longs": ["country", "division"], # essential. Maybe False.
    "lat_long_defs": '../../fauna/source-data/geo_lat_long.tsv',
    "reference": {
        "path": "metadata/zika_outgroup.gb",
        "metadata": {
            'strain': "reference", "accession": "KX369547", "date": "2013-10-25",
            'host': "human", 'country': "French Polynesia"
        },
        "use": False,
        "genes": ['CA', 'PRO', 'MP', 'ENV', 'NS1', 'NS2A',
                  'NS2B', 'NS3', 'NS4A', 'NS4B', 'NS5']
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
