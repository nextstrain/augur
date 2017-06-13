from __future__ import print_function
import os, sys
sys.path.append('..')
from base.prepare import prepare
from datetime import datetime
from base.utils import fix_names
import argparse

def collect_args():
    parser = argparse.ArgumentParser(description = "Prepare fauna FASTA for analysis")
    parser.add_argument('-s', '--serotype', choices=["all", "denv1", "denv2", "denv3", "denv4"], default="all", type=str, help="serotype (default: all)")
    return parser.parse_args()

dropped_strains = [
    'DENV1/VIETNAM/BIDV992/2006', 'DENV1/FRANCE/00475/2008', 'DENV1/VIETNAM/BIDV3990/2008', 'DENV2/HAITI/DENGUEVIRUS2HOMOSAPIENS1/2016', # Probable recombinants
    'DENV2/AUSTRALIA/QML22/2015' # Suspiciously far diverged
]
references = {
    "denv1": {"metadata": {'strain': "DENV1/NAURUISLAND/REFERENCE/1997", "accession": "NC_001477", "date": "1997-XX-XX", 'host': "NA", 'country': "Nauru"}},
    "denv2": {"metadata": {'strain': "DENV2/THAILAND/REFERENCE/1964", "accession": "NC_001474", "date": "1964-XX-XX", 'host': "NA", 'country': "Thailand"}},
    "denv3": {"metadata": {'strain': "DENV3/SRI_LANKA/REFERENCE/2000", "accession": "NC_001475", "date": "2000-XX-XX", 'host': "NA", 'country': "Sri Lanka"}},
    "denv4": {"metadata": {'strain': "DENV4/NA/REFERENCE/2003", "accession": "NC_002640", "date": "2003-XX-XX", 'host': "NA", 'country': "NA"}},
}
for key in references:
    references[key]["genes"] = ['C', 'M', 'E', 'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', '2K', 'NS4B', 'NS5']
    references[key]["use"] = False
    references[key]["path"] = "metadata/dengue_{}_outgroup.gb".format(key)
references["all"] = references["denv4"]

config = {
    "dir": "dengue",
    "file_prefix": "dengue_all",
    "input_paths": [],
    "header_fields": {0:'strain', 1:'accession', 2:'date', 3:'region', 4:'country',
                    5:'division', 6: 'location', 7: 'authors', 8: 'url'},
    "filters": (
        ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
        ("Restrict Date Range", lambda s: s.attributes['date'] >= datetime(2012,01,1).date()),
        ("Restrict Date Range", lambda s: s.attributes['date'] <= datetime(2018,01,1).date()),
        ("Sequence Length", lambda s: len(s.seq)>=5000),
        ("Bad Region", lambda s: s.attributes['region'] not in ['', '?'])
    ),
    "subsample": {
        "category": lambda x:(x.attributes['region'], x.attributes['date'].year, x.attributes['date'].month),
        "threshold": 3,
    },
    "add_urls": {
        "prefix": "https://www.ncbi.nlm.nih.gov/nuccore/%s",
        "attr": "accession"
    },
    "colors": ["region"], # essential. Maybe False.
    "color_defs": ["./colors.tsv"],
    "lat_longs": ["region"], # essential. Maybe False.
    "lat_long_defs": '../../fauna/source-data/geo_lat_long.tsv',
    "reference": {
        "path": "metadata/ebola_outgroup.gb",
        "metadata": {
            'strain': "reference", "accession": "KR075003", "date": "2014-XX-XX",
            'host': "human", 'country': "Liberia"
        },
        "include": 0,
        "genes": ['NP', 'VP35', 'VP40', 'GP', 'sGP', 'VP30', 'VP24', 'L']
    }
}


if __name__=="__main__":
    params = collect_args()

    # modify config file according to serotype
    config["input_paths"] = ["../../fauna/data/dengue_{}.fasta".format(params.serotype)]
    config["reference"] = references[params.serotype]
    config["file_prefix"] = "dengue_{}".format(params.serotype)

    runner = prepare(config)
    runner.load_references()
    runner.applyFilters()
    runner.ensure_all_segments()
    runner.subsample()
    runner.colors()
    runner.latlongs()
    runner.write_to_json()
