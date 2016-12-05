from __future__ import division, print_function
from collections import defaultdict
import sys
sys.path.append('')  # need to import from base
from base.process import process
import numpy as np
from datetime import datetime
from base.io_util import myopen

regions = ['africa', 'south_asia', 'europe', 'china', 'north_america',
           'china', 'south_america', 'japan_korea', 'oceania', 'southeast_asia']

region_groups = {'NA':'north_america',
                 'AS':['china', 'japan_korea', 'south_asia', 'southeast_asia'],
                 'OC':'oceania', 'EU':'europe'}

attribute_nesting = {'geographic location':['region', 'country', 'city'],}


def sampling_category(x):
    return (x.attributes['region'],
            x.attributes['date'].year,
            x.attributes['date'].month)


def sampling_priority(seq):
    return len(seq.seq)*0.0001 - 0.01*np.sum([seq.seq.count(nuc) for nuc in 'NRWYMKSHBVD'])


if __name__=="__main__":
    import argparse
    import matplotlib.pyplot as plt
    plt.ion()

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('-v', '--viruses_per_month', type = int, default = 10, help='number of viruses sampled per month')
    parser.add_argument('-y', '--resolution', type = str, default = '3y', help='outfile suffix')
    parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
    parser.add_argument('-d', '--download', action='store_true', default = False, help='load from database')
    parser.add_argument('-t', '--time_interval', nargs=2, default=('2012-01-01', '2016-01-01'),
                            help='time interval to sample sequences from: provide dates as YYYY-MM-DD')
    parser.add_argument('-l', '--lineage', type = str, default = 'h3n2', help='flu lineage to process')
    parser.add_argument('--load', action='store_true', help = 'recover from file')
    parser.add_argument('--no_tree', default=False, action='store_true', help = "don't build a tree")
    params = parser.parse_args()

    ppy = 12
    time_interval = [datetime.strptime(x, '%Y-%m-%d').date()  for x in params.time_interval]
    pivots = np.arange(time_interval[0].year+(time_interval[0].month-1)/12.0,
                       time_interval[1].year+time_interval[1].month/12.0, 1.0/ppy)

    # load data from all segments
    segment_names = ['pb1', 'pb2', 'pa', 'ha', 'np', 'na', 'm', 'ns']
    segments = {}
    viruses = defaultdict(list)
    for seg in segment_names:
        input_data_path = '../fauna/data/'+params.lineage+'_'+seg
        if seg=='m':
            input_data_path+='p'
        store_data_path = 'store/'+params.lineage + '_' + params.resolution + '_' + seg + '_'
        build_data_path = 'build/'+params.lineage + '_' + params.resolution + '_' + seg + '_'
        flu = process(input_data_path = input_data_path, store_data_path = store_data_path,
                       build_data_path = build_data_path, reference='flu/metadata/'+params.lineage + '_' + seg +'_outgroup.gb',
                       proteins=['SigPep', 'HA1', 'HA2'],
                       method='SLSQP', inertia=np.exp(-1.0/ppy), stiffness=2.*ppy)

        flu.load_sequences(fields={0:'strain', 2:'isolate_id', 3:'date', 4:'region',
                             5:'country', 7:"city", 12:"subtype",13:'lineage'})

        print("## loading data for segment %s, found %d number of sequences"%(seg, len(flu.seqs.all_seqs)))
        for sequence in flu.seqs.all_seqs:
            viruses[sequence].append(seg)

        segments[seg] = flu

    # determine strains that are complete
    complete_strains = filter(lambda x:len(viruses[x])==len(segment_names), viruses.keys())
    # determine filter every segment down to the sequences for which all other segments exist
    segments['ha'].seqs.filter(lambda s: s.name in complete_strains)
    segments['ha'].seqs.filter(lambda s:s.attributes['date']>=time_interval[0] and s.attributes['date']<time_interval[1])
    segments['ha'].seqs.subsample(category = sampling_category, priority = sampling_priority, threshold = params.viruses_per_month)
    strains_to_use = segments['ha'].seqs.seqs.keys()

    # align and build tree
    for seg, flu in segments.iteritems():
        flu.seqs.filter(lambda s: s.name in strains_to_use)
        flu.seqs.filter(lambda s:
            s.attributes['date']>=time_interval[0] and s.attributes['date']<time_interval[1])
        if seg!='ha':
            flu.seqs.seqs = flu.seqs.all_seqs

        flu.align()
        flu.dump()
        flu.build_tree()
        flu.annotate_tree(Tc=0.005, timetree=True, reroot='best')
        flu.tree.geo_inference('region')

        flu.dump()
        flu.export(extra_attr=[], controls=attribute_nesting)

    # determine ladder rank strains in of every tree
    ladder_ranks = defaultdict(list)
    for seg in segment_names:
        for leaf in segments[seg].tree.tree.get_terminals():
            ladder_ranks[leaf.name].append(leaf.yvalue)


    for seg in segment_names:
        for leaf in segments[seg].tree.tree.get_terminals():
            leaf.attr['ladder_ranks'] = list(np.array(ladder_ranks[leaf.name]))
            leaf.nleafs=1

    for seg in segment_names:
        for node in segments[seg].tree.tree.get_nonterminals(order='postorder'):
            node.nleafs = np.sum([x.nleafs for x in node])
            node.attr['ladder_ranks'] = list(np.sum([x.nleafs*x.attr['ladder_ranks'] for x in node], axis=0)/node.nleafs)

    # align and build tree
    for seg, flu in segments.iteritems():
        flu.dump()
        flu.export(extra_attr=[], controls=attribute_nesting)
