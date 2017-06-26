"""
This script takes fauna fasta files and produces filtered and subsetted JSONs
To be processed by augur
There are some custom command-line arguments - I envisage this will be somewhat
common for different analyses
"""
from __future__ import print_function
import os, sys, re
sys.path.append('..') # we assume (and assert) that this script is running from inside the flu folder
from base.prepare import prepare
from datetime import datetime, timedelta, date
from base.utils import fix_names
import argparse
from pprint import pprint
from pdb import set_trace
from flu_info import regions, outliers, reference_maps, reference_viruses, segments
from flu_subsampling import flu_subsampling

def collect_args():
    parser = argparse.ArgumentParser(description = "Prepare fauna FASTA for analysis")
    parser.add_argument('-l', '--lineages', choices=['h3n2', 'h1n1pdm', 'vic', 'yam'], default=['h3n2', 'h1n1pdm', 'vic', 'yam'], nargs='+', type=str, help="serotype (default: 'h3n2', 'h1n1pdm', 'vic' & 'yam')")
    parser.add_argument('-r', '--resolutions', default=['3y', '6y', '12y'], nargs='+', type = str,  help = "resolutions (default: 3y, 6y & 12y)")
    parser.add_argument('-s', '--segments', choices=segments, default=['ha'], nargs='+', type = str,  help = "segments (default: ha)")
    parser.add_argument('-v', '--viruses_per_month_seq', type = int, default = 0, help='number of viruses sampled per month (optional) (defaults: 90 (2y), 60 (3y), 24 (6y) or 12 (12y))')
    parser.add_argument('--sampling', default = 'even', type=str,
                        help='sample evenly over regions (even) (default), or prioritize one region (region name), otherwise sample randomly')
    parser.add_argument('--strains', help="a list of specific strains to prepare without filtering or subsampling")
    return parser.parse_args()

# for flu, config is a function so it is applicable for multiple lineages
def make_config(lineage, resolution, params):
    years_back = int(re.search("(\d+)", resolution).groups()[0])
    time_interval = [datetime.strptime(x, '%Y-%m-%d').date() for x in ["{:%Y-%m-%d}".format(datetime.today()), "{:%Y-%m-%d}".format(datetime.today()  - timedelta(days=365.25 * years_back))]]
    reference_cutoff = date(year = time_interval[0].year - 3, month=1, day=1)

    return {
        "dir": "flu",
        "file_prefix": "flu_{}".format(lineage),
        "segments": params.segments,
        "resolution": resolution,
        "lineage": lineage,
        "input_paths": ["../../fauna/data/{}_{}.fasta".format(lineage, segment) for segment in params.segments],
        #  0                     1   2         3          4      5     6       7       8          9                             10  11
        # >A/Galicia/RR9542/2012|flu|EPI376225|2012-02-23|europe|spain|galicia|galicia|unpassaged|instituto_de_salud_carlos_iii|47y|female
        "header_fields": {
            0:'strain',  2:'isolate_id', 3:'date',
            4:'region',  5:'country',    6:'division',
            8:'passage', 9:'lab',        10:'age',
            11:'gender'
        },
        "filters": (
            ("Time Interval", lambda s:
                (s.attributes['date']<=time_interval[0] and s.attributes['date']>=time_interval[1]) or
                (s.name in reference_viruses[lineage] and s.attributes['date']>reference_cutoff)
            ),
            ("Sequence Length", lambda s: len(s.seq)>=900),
            # what's the order of evaluation here I wonder?
            ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in outliers[lineage]]),
            ("Bad geo info", lambda s: s.attributes["country"]!= "?" and s.attributes["region"]!= "?" ),
        ),
        "subsample": flu_subsampling(params, years_back, "../../fauna/data/{}_crick_hi_cell".format(lineage)),
        "colors": ["region"],
        "color_defs": ["colors.flu.tsv"],
        "lat_longs": ["country", "region"],
        "lat_long_defs": '../../fauna/source-data/geo_lat_long.tsv',
        "references": {seg:reference_maps[lineage][seg] for seg in params.segments},
        "regions": regions,
        "time_interval": time_interval,
        "strains": params.strains
    }

if __name__=="__main__":
    params = collect_args()
    # set_trace()

    ## lots of loops to allow multiple downstream analysis
    for lineage in params.lineages:
        for resolution in params.resolutions:
            pprint("Preparing lineage {}, segments: {}, resolution: {}".format(lineage, params.segments, resolution))

            config = make_config(lineage, resolution, params)
            runner = prepare(config)
            runner.load_references()
            if params.strains is None:
                runner.applyFilters()
            runner.ensure_all_segments()
            runner.subsample()
            runner.colors()
            runner.latlongs()
            runner.write_to_json()
