from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import base.prepare
from base.prepare import prepare
from datetime import datetime
from base.utils import fix_names
import argparse

dropped_strains = [
    'KAK_428', 'NGA/2016/ISTH_0779', 'GUINEA_Z_185a', 'NIGERIA_IKEJI' # overly divergent
]
forced_strains = [
]

def collect_args():
    """Returns an Ebola-specific argument parser.
    """
    parser = base.prepare.collect_args()
    segments = ["s", "l"]
    parser.add_argument('-s', '--segments', choices=segments, default=segments, nargs='+', type = str,  help = "segments to prepare")
    parser.set_defaults(
        viruses_per_month=0
    )
    return parser

def make_config(params):
    return {
        "dir": "lassa",
        "file_prefix": "lassa",
        "title": "Real-time tracking of Lassa virus evolution",
        "maintainer": ["Bedford Lab", "http://bedford.io/team/"],
        "input_paths": [
            "../../../flora/data/lassa_s.fasta",
            "../../../flora/data/lassa_l.fasta",
        ],
        "header_fields": {0:'strain', 1:'accesion', 2: 'segment', 3:'date', 4:'region', 5: 'country', 6:'host_species', 7:'authors', 8:'title', 9:'journal', 10:'paper_url'},
        "filters": (
            ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
            ("Restrict Date Range for S segment", {
                "s": lambda s: True,
                "l": lambda s: True
            }),
            # ("Restrict Date Range", lambda s: s.attributes['date'] <= datetime(2018,01,1).date()),
            ("Sequence Length", {
                "s": lambda s: len(s.seq)>=2500,
                "l": lambda s: len(s.seq)>=5000,
            })
        ),
        "subsample": {
            "category": lambda x:(x.attributes['country'], x.attributes['date'].year),
            "threshold": params.viruses_per_month,
            "priority": lambda x:x.id in forced_strains
        },
        "colors": ["country", "host_species"],
        "color_defs": ["./colors.tsv"],
        "lat_longs": ["country"],
        "auspice_filters": ["country", "authors", "host_species"],
        "references": {
            # references are pinneo strain. Same as Kristian's Cell paper.
            # Pinneo paper: http://jvi.asm.org/content/74/15/6992.long
            # Cell paper: http://www.cell.com/cell/pdfExtended/S0092-8674(15)00897-1
            "s": {
                "path": "metadata/lassa_s.gb",
                "metadata": {
                    'strain': "Nig08_04", "accession": "GU481068", "date": "2008-XX-XX",
                    'country': "nigeria", 'segment': 'S'
                },
                "include": 1,
                "genes": ['NP', 'GPC']
            },
            "l": {
                "path": "metadata/lassa_l.gb",
                "metadata": {
                    'strain': "Pinneo-NIG-1969", "accession": "KM822127", "date": "1969-XX-XX",
                    'country': "nigeria", 'segment': 'L'
                },
                "include": 1,
                "genes": ['Z', 'L']
            },
        }
    }


if __name__=="__main__":
    parser = collect_args()
    params = parser.parse_args()
    config = make_config(params)
    config["segments"] = params.segments

    if params.viruses_per_month == 0:
        config["subsample"] = False

    runner = prepare(config)
    runner.load_references()
    runner.applyFilters()
    # runner.ensure_all_segments()
    runner.subsample()
    runner.colors()
    runner.latlongs()
    runner.write_to_json(segment_addendum=True)
