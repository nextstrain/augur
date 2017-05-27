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

regions = ['africa', 'south_asia', 'europe', 'china', 'north_america',
           'china', 'south_america', 'japan_korea', 'oceania', 'southeast_asia']

outliers = {
'h3n2':["A/Sari/388/2006", "A/SaoPaulo/36178/2015", "A/Pennsylvania/40/2010", "A/Pennsylvania/14/2010",
        "A/Pennsylvania/09/2011", "A/OSAKA/31/2005", "A/Ohio/34/2012", "A/Kenya/170/2011", "A/Kenya/168/2011",
        "A/Indiana/21/2013", "A/Indiana/13/2012", "A/Indiana/6/2013", "A/Indiana/17/2013", "A/Indiana/11/2013", "A/Indiana/08/2012", "A/Indiana/06/2013",
        "A/India/6352/2012", "A/HuNan/01/2014", "A/Helsinki/942/2013", "A/Guam/AF2771/2011", "A/Chile/8266/2003",
        "A/Busan/15453/2009", "A/Nepal/142/2011", "A/Kenya/155/2011", "A/Guam/AF2771/2011", "A/Michigan/82/2016",
        "A/Ohio/27/2016", "A/Ohio/28/2016", "A/Michigan/83/2016", "A/Michigan/84/2016", "A/Jiangsu-Tianning/1707/2013",
        "A/HuNan/1/2014", "A/Iran/227/2014", "A/Iran/234/2014", "A/Iran/140/2014", "A/Jiangsu-Chongchuan/1830/2014",
        "A/Chile/8266/2003", "A/Louisiana/4/2003", "A/Lousiana/4/2003", "A/OSAKA/31/2005",
        "A/Sari/388/2006", "A/HongKong/HK1/2008", "A/HongKong/HK1MA21-1/2008","A/HongKong/HK1MA21-4/2008", "A/HongKong/HK1MA21-2/2008",
        "A/HongKong/HK1MA21-3/2008", "A/HongKong/HK2/2008", "A/HongKong/HK2MA21-1/2008",
        "A/HongKong/HK2MA21-2/2008", "A/HongKong/HK2MA21-3/2008", "A/HongKong/HK4/2008",
        "A/HongKong/HK5/2008", "A/HongKong/HK5MA21-1/2008", "A/HongKong/HK5MA21-3/2008",
        "A/HongKong/HK6/2008", "A/HongKong/HK6MA21-2/2008", "A/HongKong/HK6MA21-3/2008",
        "A/HongKong/HK7/2008", "A/HongKong/HK8/2008", "A/HongKong/HK8MA21-1/2008",
        "A/HongKong/HK8MA21-2/2008", "A/HongKong/HK8MA21-3/2008", "A/HongKong/HK8MA21-4/2008",
        "A/HongKong/HK9/2008", "A/HongKong/HK9MA21-1/2008", "A/HongKong/HK9MA21-2/2008",
        "A/HongKong/HK9MA21-3/2008", "A/HongKong/HK10/2008", "A/HongKong/HK10MA21-1/2008",
        "A/HongKong/HK10MA21-2/2008", "A/HongKong/HK10MA21-3/2008", "A/HongKong/HK10MA21-4/2008",
        "A/HongKong/HK11MA21-1/2008", "A/HongKong/HK11MA21-3/2008", "A/HongKong/HK11MA21-4/2008",
        "A/HongKong/HK12/2008", "A/HongKong/HK12MA21-2/2008", "A/HongKong/HKMA12/2008",
        "A/HongKong/HKMA12A/2008", "A/HongKong/HKMA12B/2008", "A/HongKong/HKMA12D/2008",
        "A/HongKong/HKMA12E/2008", "A/HongKong/HKMA20B/2008", "A/HongKong/HKMA20E/2008", "A/Kansas/13/2009",
        "A/Busan/15453/2009", "A/Pennsylvania/14/2010", "A/Pennsylvania/40/2010", "A/Guam/AF2771/2011",
        "A/Indiana/8/2011", "A/Kenya/155/2011", "A/Kenya/168/2011", "A/Kenya/170/2011", "A/Nepal/142/2011",
        "A/Pennsylvania/09/2011", "A/Pennsylvania/9/2011", "A/Quebec/167936/2011", "A/Quebec/170658/2011",
        "A/India/6352/2012", "A/Indiana/08/2012", "A/Indiana/13/2012", "A/Ohio/34/2012",
        "A/Helsinki/942/2013", "A/Indiana/06/2013", "A/Indiana/11/2013", "A/Indiana/21/2013",
        "A/Jiangsu-Tianning/1707/2013", "A/HuNan/01/2014", "A/Jiangsu-Chongchuan/1830/2014",
        "A/Jiangsu-Chongchuan/12179/2014", "A/Ohio/2/2014", "A/Ohio/4319/2014", "A/SaoPaulo/3-34891/2014",
        "A/Wisconsin/24/2014", "A/NewJersey/53/2015", "A/SaoPaulo/36178/2015", "A/SaoPaulo/61282/2015",
        "A/SaoPaulo/65194/2015", "A/Michigan/39/2015", "A/Sydney/53/2015", "A/Michigan/82/2016",
        "A/Michigan/83/2016", "A/Michigan/84/2016", "A/Michigan/87/2016", "A/Michigan/89/2016",
        "A/Michigan/90/2016", "A/Michigan/91/2016", "A/Michigan/93/2016", "A/Michigan/94/2016",
        "A/Michigan/95/2016", "A/Michigan/96/2016", "A/Ohio/27/2016", "A/Ohio/28/2016", "A/Ohio/32/2016",
        "A/Ohio/33/2016", "A/Ohio/35/2016", "A/Zhejiang-Wuxin/1300/2016", "A/Nanjing/1/2010", "A/StPetersburg/5/2009",
        "A/Cambodia/NHRCC00001/2009","A/Cambodia/NHRCC00002/2009", "A/Cambodia/NHRCC00003/2009",
        "A/Cambodia/NHRCC00006/2009", "A/Iran/407/2014"],
'h1n1pdm': ["A/Kenya/264/2012", "A/Iowa/39/2015", "A/Asturias/RR6898/2010", "A/Wisconsin/28/2011",
            "A/Brest/1161/2014", "A/Tomsk/273-MA1/2010", "A/Minnesota/46/2015", "A/Poland/16/2013",
            "A/Hungary/02/2013", "A/Hungary/16/2013", "A/California/07/2009NYMC-X18113/198",
            "A/Christchurch/16/2010NIB-74xp13/202", "A/Bari/166/2016", "A/Bari/167/2016", "A/Dakar/3/2014",
            "A/Arkansas/15/2013", "A/Wisconsin/87/2005", 'A/Norway/1929/2014', 'A/Ohio/9/2015'],
'vic':["B/Bangkok/SI17/2012", "B/Bangkok/SI58/2012", "B/Kol/2024/2008", "B/Kolkata/2024/2008"],
"yam":[]
}

reference_maps = {
    "h3n2": {
        "ha": {
            "path": "metadata/h3n2_ha_outgroup.gb",
            "metadata": {
                'strain': "A/Beijing/32/1992",
                'isolate_id': "CY113677",
                'date': "1992-XX-XX",
                'region': "china",
                'country': "China",
                "city": "Beijing",
                "passage": "unknown",
                'lab': "unknown",
                'age': "unknown",
                'gender': "unknown"
            },
            "use": True,
            "genes": ["HA1", "HA2"]
        }
    }
}

input_paths = {
    "h3n2": {
        "ha": "../../fauna/data/h3n2.fasta",
    },
}

reference_viruses = {
    'h3n2': ['A/Wisconsin/67/2005', 'A/Brisbane/10/2007',  'A/Perth/16/2009', 'A/Victoria/361/2011','A/Texas/50/2012', 'A/Switzerland/9715293/2013', 'A/HongKong/4801/2014', 'A/Alaska/232/2015'],
    'h1n1pdm':[],
    'vic':[],
    'yam':[]
}

# for flu, config is a function so it is applicable for multiple lineages
def make_config (lineage, segments, time_interval, reference_cutoff):
    return {
        "dir": "flu", # the current directory. You mush be inside this to run the script.
        "file_prefix": "flu_{}".format(lineage),
        "segments": segments,
        "input_format": "fasta",
        "input_paths": [input_paths[lineage][segment] for segment in segments],
        "output_folder": "prepared",
        # note that "strain" is essential and "date" is special
        "header_fields": {
            0:'strain',
            2:'isolate_id',
            3:'date',
            4:'region',
            5:'country',
            7:"city",
            8:"passage",
            9:'lab',
            10:'age',
            11:'gender'
        },
        # can provide multiple date formats here - FIFO
        "date_format": ["%Y-%m-%d"],
        "require_dates": True,
        # require that all selected isolates have all the genomes
        "ensure_all_segments": False,
        # see docs for help with filters - nested tuples abound
        "filters": (
            ("Time Interval", lambda s:
                (s.attributes['date']>=time_interval[0] and s.attributes['date']<time_interval[1]) or
                (s.name in reference_viruses[lineage] and s.attributes['date']>reference_cutoff)
            ),
            ("Sequence Length", lambda s: len(s.seq)>=900),
            # what's the order of evaluation here I wonder?
            ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in outliers[lineage]]),
        ),
        # see the docs for this too! if you don't want to subsample, set it to False
        # "subsample": False,
        "subsample": {
            "category": None,
            "priority": None,
            "threshold": 2,
            # somehow pass regions in here i guess

        },

        # see the docs for what's going on with colours (sic) & lat/longs
        "colors": ["country", "region", "city"], # essential. Maybe False.
        "color_defs": ["colors.flu.tsv"],
        "lat_longs": ["country", "division"], # essential. Maybe False.
        "lat_long_defs": '../../fauna/source-data/geo_lat_long.tsv',

        "references": {seg:reference_maps[lineage][seg] for seg in segments},

    }



if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "Prepare a JSON for processing")
    parser.add_argument('--bin', type = str, default = "python")
    parser.add_argument('--builds', default=['live', 'cdc'], nargs='+', type = str, help = "builds to include")
    parser.add_argument('--lineages', default=['h3n2', 'h1n1pdm', 'vic', 'yam'], nargs='+', type = str,  help = "lineages to include")
    parser.add_argument('--resolutions', default=['2y', '3y', '6y', '12y'], nargs='+', type = str,  help = "resolutions to include")
    parser.add_argument('--assays', default=['hi', 'fra'], nargs='+', type = str, help = "assays to include (CDC build only)")
    parser.add_argument('--passages', default=['egg', 'cell'], nargs='+', type = str, help = "passages to include (CDC build only)")
    parser.add_argument('--live_auspice_path', default = '../nextflu/auspice/data/')
    parser.add_argument('--cdc_auspice_path', default = '../nextflu-cdc/auspice/data/')
    params = parser.parse_args()

    # for testing:
    params.resolutions = ["3yr"]
    params.lineages = ["h3n2"]

    if 'live' in params.builds:
        for lineage in params.lineages:
            for resolution in params.resolutions:
                segments = ["ha"]

                print ('\n--------------------\n')
                print ('Processing lineage {} with resolution {} for all HI data, for segments {}'.format(lineage, resolution, ", ".join(segments)))

                # time_interval can be defined in seasonal_flu.py, but
                # run.py uses these defaults:
                years_back = int(re.search("(\d+)", resolution).groups()[0])
                today_str = "{:%Y-%m-%d}".format(datetime.today())
                date_str = "{:%Y-%m-%d}".format(datetime.today() - timedelta(days=365.25 * years_back))
                time_interval = [datetime.strptime(x, '%Y-%m-%d').date() for x in [date_str, today_str]]
                reference_cutoff = date(year = time_interval[0].year - 3, month=1, day=1)
                config = make_config(lineage, segments, time_interval, reference_cutoff)

                pprint("SUBSAMPLING STILL TO DO")
                # pprint(config)

                runner = prepare(config)
                runner.load_references()
                runner.applyFilters()
                runner.ensure_all_segments()
                runner.subsample()
                runner.colors()
                runner.latlongs()
                runner.write_to_json()
