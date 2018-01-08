from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
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

    segments = ['pb2', 'pb1', 'pa', 'ha', 'np', 'na', 'mp', 'ns']
    parser.add_argument('-l', '--lineage', choices=['h7n9'], default='h7n9', type=str, help="subtype")
    parser.add_argument('-s', '--segments', choices=segments, default=segments, nargs='+', type = str,  help = "segments to prepare")
    parser.set_defaults(
        viruses_per_month=0
    )

    return parser

dropped_strains = [
    "A/Chicken/Netherlands/16007311-037041/2016", "A/Chicken/Netherlands/1600",
    "A/duck/Zhejiang/LS02/2014", "A/northernshoveler/California/AH0076755/2016",
    "A/chicken/Tennessee/17-008279-6/2017", "A/duck/Alabama/17-008643-2/2017",
    "A/Americangreen-wingedteal/Idaho/AH0006121/2015",
    "A/cinnamonteal/California/AH0079011/2016", "A/chicken/Tennessee/17-007147-7/2017",
    "A/chicken/Alabama/17-008899-11/2017", "A/chicken/Alabama/17-008899-10/2017",
    "A/Duck/Bangladesh/2698", "A/guineafowl/Alabama/17-008646-1/2017",
    "A/chicken/Tennessee/17-007429-7/2017", "A/Duck/Bangladesh/2699",
    "A/nothernshoveler/California/AH0079371/2016", "A/duck/Shanghai/SD016/2015",
    "A/duck/Zhejiang/S4488/2014", "A/duck/Shanghai/SD015/2015",
        # not part of epi clade
    "A/British Columbia/1/2015", "A/BritishColumbia/1/2015",
    "A/blue-wingedteal/Louisiana/UGAI15-1692/2015",
    "A/blue-wingedteal/Louisiana/UGAI15-1367/2015",
        # travel case. throws off map
    "A/Duck/Bangladesh/2704", # clock is off
    "_A/Malaysia/228/2014_H7N9__", # travel case
]

config = {
    "dir": "avian", # the current directory. You must be inside this to run the script.
    "title": "Real-time tracking of influenza A/H7N9 evolution",
    "maintainer": ["James Hadfield", "http://bedford.io/team/james-hadfield/"],
    "input_paths": [
        "../../../fauna/data/h7n9_pb2.fasta",
        "../../../fauna/data/h7n9_pb1.fasta",
        "../../../fauna/data/h7n9_pa.fasta",
        "../../../fauna/data/h7n9_ha.fasta",
        "../../../fauna/data/h7n9_np.fasta",
        "../../../fauna/data/h7n9_na.fasta",
        "../../../fauna/data/h7n9_mp.fasta",
        "../../../fauna/data/h7n9_ns.fasta"
    ],
    "header_fields": {
        #  0                   1    2         3          4     5     6     7       8       9   10                                    11
        # >A/Jiangsu/9389/2014|h7n9|EPI628192|2014-01-07|human|china|china|jiangsu|jiangsu|egg|who_chinese_national_influenza_center|2017-03-10
        0:'strain', 2:'accession', 3:'date', 4:'host', 6:'country',
        7: 'division', 10: 'authors', 11: 'fauna_date'
    },
    "traits_are_dates": ["fauna_date"],
    "filters": (
        ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
        ("Prior to Epidemic", lambda s: s.attributes['date'] >= datetime(2013,1,1).date()),
        ("Missing Month Data", lambda s: "-XX-XX" not in s.attributes['raw_date']),
        ("Exclude bad host", lambda s: s.attributes["host"] not in ["laboratoryderived", "watersample"]),
        ("Not part of epi clade", lambda s: s.attributes["division"] not in ["alabama", "tennessee"]),
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
    "colors": ["division", "host"], # essential. Maybe False.
    "color_defs": ["./colors.avian.tsv"],
    "auspice_filters": ["division", "host"],
    "lat_longs": ["division"], # essential. Maybe False.
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
        segment_addendum = False
    else:
        config["file_prefix"] = "avian_h7n9"
        segment_addendum = True

    config["segments"] = params.segments

    runner = prepare(config)
    runner.load_references()
    runner.applyFilters()
    #runner.ensure_all_segments()
    runner.subsample()
    runner.colors()
    runner.latlongs()
    runner.write_to_json(segment_addendum=segment_addendum)
