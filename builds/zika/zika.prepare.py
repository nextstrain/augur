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
    parser.set_defaults(
        viruses_per_month=20,
        file_prefix="zika"
    )
    return parser

dropped_strains = [
    "CX17", # large indels in alignment
    "ZF36_36S", # possible contamination
    "Dominican_Republic/2016/PD2", "GD01", "GDZ16001", "VEN/UF_2/2016", "ZZ_1", # true strains, but duplicates of other strains in dataset
    "Bahia04", "JAM/2016/WI_JM6", "Bahia11", "Bahia12", "DOM/2016/MA_WGS16_009", "VE_Ganxian",
        "BRA/2016/FC_DQ60D1", "CX5", "Thailand/1605aTw", # excessive terminal branch length
    "VR10599/Pavia/2016", "34997/Pavia/2016", # exports
    "THA/PLCal_ZV/2013", "SK403/13AS", "SV0010/15", "SK364/13AS", # clock is off
    "MEX/42_018_F_V3_O/2016", "MEX/2016/mex50", "GTM/2016/Guatemala_2411", # clock is off
    "AFMC_U", "AFMC_S", "JMB_185", "VIE/Bra/2016", "PHL/2012/CPC_0740", # too basal
    "Boracay/16423", "SI_BKK01", "DMSc05684_16", "ZK_YN001", # too basal
    "COL/PRV_00027/2015", "COL/PRV_00028/2015", "COL/PAN_00029/2015" "COL/PAN_00030/2015",
    "COL/FLR_00034/2015", "COL/FLR_00035/2015", "COL/FLR_00036/2015", "COL/FLR_00038/2015",
    "COL/FLR_00040/2015", "COL/FLR_00041/2015", "COL/FLR_00042/2015" # almost certainly misdated, other Colombian genomes from same study are Dec 2015 and fit the clock
]

config = {
    "dir": "zika",
    "file_prefix": "zika",
    "title": "Real-time tracking of Zika virus evolution",
    "maintainer": ["Trevor Bedford", "http://bedford.io/team/trevor-bedford/"],
    "input_paths": ["../../../fauna/data/zika.fasta"],
    "header_fields": {0:'strain', 2:'accession', 3:'date', 4:'region', 5:'country',
                    6:'division', 8:'db', 10:'authors', 11:'url', 12:'title',
                    13: 'journal', 14: 'paper_url'},
    "filters": (
        ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
        ("Restrict Date Range", lambda s: s.attributes['date'] >= datetime(2012,01,1).date()),
        ("Restrict Date Range", lambda s: s.attributes['date'] <= datetime(2018,01,1).date()),
        ("Sequence Length", lambda s: len(s.seq)>=5000),
    ),
    "subsample": {
        "category": lambda x:(x.attributes['date'].year, x.attributes['date'].month, x.attributes['country']),
    },
    "colors": ["authors", "region", "country"],
    "color_defs": ["config/colors.tsv"],
    "lat_longs": ["country", "region"],
    "auspice_filters": ["authors", "region", "country"],
    "reference": {
        "path": "config/zika_outgroup.gb",
        "metadata": {
            'strain': "PF13/251013_18", "accession": "KX369547", "date": "2013-10-25",
            'host': "human", 'country': "french_polynesia"
        },
        "include": 2,
        "genes": ['CA', 'PRO', 'MP', 'ENV', 'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', '2K', 'NS4B', 'NS5']
    }
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

    runner = prepare(config)
    runner.load_references()
    runner.applyFilters()
    runner.ensure_all_segments()
    runner.subsample()
    runner.colors()
    runner.latlongs()
    runner.write_to_json()
