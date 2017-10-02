from __future__ import print_function
import os, sys, re
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import base.prepare
from base.prepare import prepare
from base.titer_model import TiterModel
from datetime import datetime, timedelta, date
from base.utils import fix_names
from pprint import pprint
from pdb import set_trace
from flu_info import regions, outliers, reference_maps, reference_viruses, segments
from flu_subsampling import flu_subsampling

import logging

def parse_args():
    """Returns a seasonal flu-specific argument parser.
    """
    parser = base.prepare.collect_args()

    parser.add_argument('-l', '--lineage', choices=['h3n2', 'h1n1pdm', 'vic', 'yam'], default='h3n2', type=str, help="lineage (default: h3n2)")
    parser.add_argument('-r', '--resolution', default=['3y'], nargs='+', type = str,  help = "resolutions (default: 3y)")
    parser.add_argument('-s', '--segment', choices=segments, default=['ha'], nargs='+', type = str,  help = "segment (default: ha)")
    parser.add_argument('--sampling', default = 'even', type=str,
                        help='sample evenly over regions (even) (default), or prioritize one region (region name), otherwise sample randomly')
    parser.add_argument('--time_interval', nargs=2, help="explicit time interval to use -- overrides resolutions"
                                                                     " expects YYYY-MM-DD YYYY-MM-DD")
    parser.add_argument('--strains', help="a text file containing a list of strains (one per line) to prepare without filtering or subsampling")
    parser.add_argument('--titers', help="tab-delimited file of titer strains and values from fauna (e.g., h3n2_hi_titers.tsv)")
    parser.add_argument('--complete_frequencies', action='store_true', help="compute mutation frequences from entire dataset")

    parser.set_defaults(
        viruses_per_month=0
    )

    return parser.parse_args()

def make_title(lineage, resolution):
    flu_type = "B" if lineage in ["vic", "yam"] else "A"
    prettyLineage = ""
    if flu_type == "A":
        prettyLineage = lineage.upper()
    if flu_type == "B":
        prettyLineage = lineage.capitalize()
    return "Real-time tracking of influenza {}/{} evolution".format(flu_type, prettyLineage)

# for flu, config is a function so it is applicable for multiple lineages
def make_config(lineage, resolution, params):
    years_back = int(re.search("(\d+)", resolution).groups()[0])
    if params.time_interval:
        time_interval = sorted([datetime.strptime(x, '%Y-%m-%d').date() for x in params.time_interval], reverse=True)
    else:
        time_interval = [datetime.today().date(), (datetime.today()  - timedelta(days=365.25 * years_back)).date()]
    reference_cutoff = date(year = time_interval[1].year - 3, month=1, day=1)
    fixed_outliers = [fix_names(x) for x in outliers[lineage]]
    fixed_references = [fix_names(x) for x in reference_viruses[lineage]]

    if params.titers is not None:
        titer_values, strains, sources = TiterModel.load_from_file(params.titers)
    else:
        titer_values = None

    if params.sequences is not None:
        input_paths = params.sequences
    else:
        input_paths = ["../../../fauna/data/{}_{}.fasta".format(lineage, segment) for segment in params.segment]

    if params.file_prefix:
        file_prefix = params.file_prefix
    else:
        file_prefix = "flu_{}_{}_{}".format(lineage, params.segment[0], resolution)

    return {
        "dir": "flu",
        "file_prefix": file_prefix,
        "title": make_title(lineage, resolution),
        "maintainer": ["Trevor Bedford", "http://bedford.io/team/trevor-bedford/"],
        "auspice_filters": ["region", "country"],
        "segments": params.segment,
        "lineage": lineage,
        "input_paths": input_paths,
        #  0                     1   2         3          4      5     6       7       8          9                             10  11
        # >A/Galicia/RR9542/2012|flu|EPI376225|2012-02-23|europe|spain|galicia|galicia|unpassaged|instituto_de_salud_carlos_iii|47y|female
        "header_fields": {
            0:'strain',  2:'isolate_id', 3:'date',
            4:'region',  5:'country',    6:'division',
            8:'passage', 9:'authors', 10:'age',
            11:'gender'
        },
        "filters": (
            ("Time Interval", lambda s:
                (s.attributes['date']<=time_interval[0] and s.attributes['date']>=time_interval[1]) or
                (s.name in fixed_references and s.attributes['date']>reference_cutoff)
            ),
            ("Sequence Length", lambda s: len(s.seq)>=900),
            # what's the order of evaluation here I wonder?
            ("Dropped Strains", lambda s: s.id not in fixed_outliers),
            ("Bad geo info", lambda s: s.attributes["country"]!= "?" and s.attributes["region"]!= "?" ),
        ),
        "subsample": flu_subsampling(params, years_back, titer_values),
        "colors": ["region"],
        "color_defs": ["colors.flu.tsv"],
        "lat_longs": ["country", "region"],
        "lat_long_defs": '../../../fauna/source-data/geo_lat_long.tsv',
        "references": {seg:reference_maps[lineage][seg] for seg in params.segment},
        "regions": regions,
        "time_interval": time_interval,
        "strains": params.strains,
        "titers": titer_values
    }

if __name__=="__main__":
    params = parse_args()
    # set_trace()
    if params.verbose:
        # Remove existing loggers.
        for handler in logging.root.handlers:
            logging.root.removeHandler(handler)

        # Set the new default logging handler configuration.
        logging.basicConfig(
            format='[%(levelname)s] %(asctime)s - %(name)s - %(message)s',
            level=logging.DEBUG,
            stream=sys.stderr
        )
        logger = logging.getLogger(__name__)
        logger.debug("Verbose reporting enabled")

    ## lots of loops to allow multiple downstream analysis
    for resolution in params.resolution:
        pprint("Preparing lineage {}, segments: {}, resolution: {}".format(params.lineage, params.segment, resolution))

        config = make_config(params.lineage, resolution, params)
        runner = prepare(config)
        runner.load_references()
        runner.applyFilters()
        runner.ensure_all_segments()
        if not params.complete_frequencies: # if complete_frequencies, do subsampling later on in flu.process
            runner.subsample()
        taxa_to_include = list(runner.segments[params.segment[0]].get_subsampled_names(config))
        runner.segments[params.segment[0]].extras['leaves'] = taxa_to_include
        runner.colors()
        runner.latlongs()
        runner.write_to_json()
        if params.time_interval:
            break
