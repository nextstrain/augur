from __future__ import division, print_function
from collections import defaultdict
import sys
sys.path.insert(0,'.')  # need to import from base
print(sys.argv, sys.path)
from base.process import process
import numpy as np
from datetime import datetime, timedelta
from base.io_util import myopen

from seasonal_flu import *


if __name__ == '__main__':
    import argparse
    import matplotlib.pyplot as plt
    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-y', '--years_back', type = str, default = 3, help='number of years back to sample')
    parser.add_argument('--resolution', type = str, help ="outfile suffix, can determine -v and -y")
    parser.add_argument('-v', '--viruses_per_month_seq', type = int, default = 10, help='number of viruses sampled per month')
    parser.add_argument('-w', '--viruses_per_month_tree', type = int, default = 10, help='number of viruses sampled per month')
    parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
    parser.add_argument('-d', '--download', action='store_true', default = False, help='load from database')
    parser.add_argument('-t', '--time_interval', nargs=2, help='specify time interval rather than use --years_back')
    parser.add_argument('-l', '--lineage', type = str, default = 'h3n2', help='flu lineage to process')
    parser.add_argument('--new_auspice', default = False, action="store_true", help='file name for new augur')
    parser.add_argument('--confidence', default = False, action="store_true", help='evaluate confidence intervals of internal node timing')
    parser.add_argument('--sampling', default = 'even', type=str,
                        help='sample evenly (even), or prioritize one region (region), otherwise sample randomly')
    parser.add_argument('--HI', default = 'hi', type=str, help='specify titer data to use')
    parser.add_argument('--load', action='store_true', help = 'recover from file')
    parser.add_argument('--no_tree', default=False, action='store_true', help = "don't build a tree")
    params = parser.parse_args()

    # default values for --viruses_per_month and --years_back from resolution
    if params.resolution == "2y":
        params.viruses_per_month_tree = 80
        params.viruses_per_month_seq = 150
        params.years_back = 2
    elif params.resolution == "3y":
        params.viruses_per_month_tree = 50
        params.viruses_per_month_seq = 100
        params.years_back = 3
    elif params.resolution == "6y":
        params.viruses_per_month_tree = 30
        params.viruses_per_month_seq = 100
        params.years_back = 6
    elif params.resolution == "12y":
        params.viruses_per_month_tree = 20
        params.viruses_per_month_seq = 100
        params.years_back = 12

    # construct time_interval from years_back
    if not params.time_interval:
        today_str = "{:%Y-%m-%d}".format(datetime.today())
        date_str = "{:%Y-%m-%d}".format(datetime.today() - timedelta(days=365.25 * params.years_back))
        params.time_interval = [date_str, today_str]

    if params.new_auspice:
        fname_prefix = "flu_"+params.lineage+'_na'
    else:
        fname_prefix = params.lineage+'_na'

    if params.sampling!="even":
        fname_prefix+='_'+params.sampling

    if params.HI!="hi":
        fname_prefix+='_'+params.HI

    input_data_path = '../fauna/data/'+params.lineage+'_na'
    if params.resolution:
        store_data_path = 'store/'+ fname_prefix + '_' + params.resolution +'_'
        build_data_path = 'build/'+ fname_prefix + '_' + params.resolution +'_'
    else:
        store_data_path = 'store/'+ fname_prefix + '_'
        build_data_path = 'build/'+ fname_prefix + '_'

    ppy = 12
    flu = flu_process(input_data_path = input_data_path, store_data_path = store_data_path,
                   build_data_path = build_data_path, reference='flu/metadata/'+params.lineage+'_na_outgroup.gb',
                   proteins=['NA'], titer=params.HI, segment='na',
                   method='SLSQP', inertia=np.exp(-1.0/ppy), stiffness=2.*ppy)


    if params.load:
        flu.load()
        flu.export()
    else:
        flu.load_sequences(fields={0:'strain', 2:'isolate_id', 3:'date', 4:'region',
                             5:'country', 7:"city", 8:"passage",9:'lab', 10:'age', 11:'gender'})

        time_interval = [datetime.strptime(x, '%Y-%m-%d').date() for x in params.time_interval]

        pivots = np.arange(time_interval[0].year+(time_interval[0].month-1)/12.0,
                           time_interval[1].year+time_interval[1].month/12.0, 1.0/ppy)
        flu.seqs.filter(lambda s:
            s.attributes['date']>=time_interval[0] and s.attributes['date']<time_interval[1])
        flu.seqs.filter(lambda s: len(s.seq)>=900)
        flu.seqs.filter(lambda s: s.name not in outliers[params.lineage])

        if params.sampling=='even':
            flu.subsample(params.viruses_per_month_seq, all_regions=False)
        elif params.sampling in regions:
            flu.subsample_priority_region(params.viruses_per_month_seq,
                                          priority_region=params.sampling, fraction=0.5)
        else:
            flu.subsample(params.viruses_per_month_seq, all_regions=True)

        flu.align()
        flu.dump()
        # first estimate frequencies globally, then region specific
        flu.estimate_mutation_frequencies(region="global", pivots=pivots, min_freq=.10)
        # for region in region_groups.iteritems():
        #     flu.estimate_mutation_frequencies(region=region)

        if not params.no_tree:
            if params.sampling=='even':
                flu.subsample(params.viruses_per_month_seq, all_regions=False, repeated=True)
            elif params.sampling in regions:
                flu.subsample_priority_region(params.viruses_per_month_tree, priority_region=params.sampling,
                                              fraction=0.5, repeated=True)
            else:
                flu.subsample(params.viruses_per_month_seq, all_regions=True, repeated=True)

            flu.align()
            flu.build_tree()
            flu.clock_filter(n_iqd=3, plot=False, remove_deep_splits=True)
            flu.annotate_tree(Tc=0.03, timetree=True, reroot='best', confidence=params.confidence)
            for geo in geo_attributes:
                flu.tree.geo_inference(geo)

            flu.estimate_tree_frequencies()
            for region in regions:
                flu.estimate_tree_frequencies(region=region)
            flu.dump()

            flu.HI_model(criterium = lambda x:len(x.aa_mutations['NA'])>0)
            flu.dump()
            flu.export(extra_attr=['serum'], controls=attribute_nesting, geo_attributes=geo_attributes,
                       color_options=color_options, panels=panels)
            flu.HI_export()
