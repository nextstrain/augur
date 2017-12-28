from __future__ import print_function
import os, sys, glob
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import base.process
from base.fitness_model import process_predictor_args
from base.process import process
from base.utils import fix_names
from flu_titers import HI_model, HI_export, H3N2_scores
from flu_info import clade_designations
import argparse
import numpy as np
from pprint import pprint
from pdb import set_trace
from Bio.Align import MultipleSeqAlignment


def parse_args():
    """Returns a seasonal flu-specific argument parser.
    """
    parser = base.process.collect_args()

    parser.add_argument('--no_mut_freqs', default=False, action='store_true', help="skip mutation frequencies")
    parser.add_argument('--no_tree_freqs', default=False, action='store_true', help="skip tree (clade) frequencies")
    parser.add_argument('--pivot_spacing', type=float, default=1.0, help="month per pivot")
    parser.add_argument('--titers_export', default=False, action='store_true', help="export titers.json file")
    parser.add_argument('--annotate_fitness', default=False, action='store_true', help="run fitness prediction model and annotate fitnesses on tree nodes")
    parser.add_argument('--predictors', default=['cTiter'], nargs='+', help="attributes to use as fitness model predictors")
    parser.add_argument('--predictors_params', type=float, nargs='+', help="precalculated fitness model parameters for each of the given predictors")
    parser.add_argument('--predictors_sds', type=float, nargs='+', help="precalculated global standard deviations for each of the given predictors")
    parser.add_argument('--epitope_mask_version', default="wolf", help="name of the epitope mask that defines epitope mutations")

    parser.set_defaults(
        json="prepared/flu.json"
    )

    return parser.parse_args()

def make_config (prepared_json, args):
    # Create a fitness model parameters data structure for the given arguments.
    predictors = process_predictor_args(
        args.predictors,
        args.predictors_params,
        args.predictors_sds
    )

    return {
        "dir": "flu",
        "in": prepared_json,
        "geo_inference": ['region'],
        "auspice": { ## settings for auspice JSON export
            "panels": ['tree', 'entropy'],
            "extra_attr": ['serum'],
            "color_options": {
                "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
            },
            "controls": {'authors':['authors']},
            "defaults": {'colorBy': 'cTiter',
                'geoResolution': 'region',
                'distanceMeasure': 'div',
                'mapTriplicate': True},
            "titers_export": args.titers_export
        },
        "titers": {
            "criterium": lambda x: len(x.aa_mutations['HA1']+x.aa_mutations['HA2'])>0,
            "lam_avi":2.0,
            "lam_pot":0.3,
            "lam_drop":2.0
        },
        "build_tree": not args.no_tree,
        "estimate_mutation_frequencies": not args.no_mut_freqs,
        "estimate_tree_frequencies": not args.no_tree_freqs,
        "epitope_mask": "metadata/h3n2_epitope_masks.tsv",
        "epitope_mask_version": args.epitope_mask_version,
        "annotate_fitness": args.annotate_fitness,
        "predictors": predictors,
        "clean": args.clean,
        "pivot_spacing": args.pivot_spacing/12.0,
        "timetree_options": {
            "Tc": 0.03,
            # "confidence":True,
            # "use_marginal":True
        },
        "newick_tree_options":{
            "raxml": not args.no_raxml
        }
    }


def rising_mutations(freqs, counts, genes, region='NA', dn=5, offset=0, baseline = 0.01, fname='tmp.txt'):
    '''
    safe a file containing all mutations and summary of their recent frequency trajectories.
    mutations are sorted by their log derivative over the past dn month
    '''
    dx = {}
    npoints = counts[region].shape[0]
    ind = np.arange(npoints)[-(dn+offset):npoints-offset]
    for gene in genes:
        for mut,f  in freqs[(region, gene)].iteritems():
            c = np.sum(counts[region][ind]*f[ind])
            tmp_x = f[ind].mean()
            tmp_dx = f[npoints-1-offset] - f[-dn-offset]
            dx[(region, gene, mut[0], mut[1])] = (tmp_x, tmp_dx, tmp_dx/(tmp_x+baseline), c)

    with open(fname, 'w') as ofile:
        ofile.write("#Frequency change over the last %d month in region %s\n"%(dn, region))
        print("#Frequency change over the last %d month in region %s"%(dn, region))
        for k,v in sorted(dx.items(), key=lambda x:x[1][2], reverse=True):
            ofile.write("%s:%d%s\t%1.3f\t%1.3f\t%1.3f\t%1.1f\n"%(k[1], k[2]+1, k[3], v[0], v[1], v[2], v[3]))
            print("%s:%d%s: x=%1.3f, dx=%1.3f, dx/x=%1.3f, c=%1.1f"%(k[1], k[2]+1, k[3], v[0], v[1], v[2], v[3]))
        ofile.write('\n')
    return dx


def recurring_mutations(tree, fname_by_position='tmp.txt', fname_by_mutation='tmp.txt'):
    '''
    count the number of times that each position has mutated on the tree and save to file.
    in additition, count the number of mutations that resulted in a specific substitution
    '''
    from collections import defaultdict
    by_mutation = defaultdict(int)
    by_position = defaultdict(int)
    for n in tree.find_clades():
        for gene in n.aa_mutations:
            for mut in n.aa_mutations[gene]:
                by_position[(gene, mut[1])]+=1
                by_mutation[(gene, mut)]+=1

    with open(fname_by_position, 'w') as ofile:
        ofile.write("#Number of independent mutation events in the tree observed at a specific position\n")
        for mut,val in sorted(by_position.items(), key=lambda x:x[1], reverse=True):
            ofile.write("%s:%d\t%d\n"%(mut[0], mut[1]+1, val))

    with open(fname_by_mutation, 'w') as ofile:
        ofile.write("#Number of sepcific independent mutation events in the tree observed\n")
        for mut,val in sorted(by_mutation.items(), key=lambda x:x[1], reverse=True):
            ofile.write("%s:%s%d%s\t%d\n"%(mut[0], mut[1][0], mut[1][1]+1, mut[1][2], val))

    return by_mutation, by_position


def freq_auto_corr(freq1, freq2, min_dfreq=0.2):
    '''
    calculate the autocorrelation function of two sets of frequencies, say in Asia and Oceania,
    to see whether one precedes the other on average. EXPERIMENTAL
    '''
    dt = 5
    corr = np.zeros(2*dt+1)
    for mut in freq1:
        if mut in freq2:
            f1 = freq1[mut]
            if f1.max() - f1.min()>min_dfreq and f1[-1]>0.8 and f1[0]<0.2:
                f2 = freq2[mut]
                corr[dt]+= np.mean((f1-f1.mean())*(f2-f2.mean()))
                for i in range(1,dt+1):
                    corr[dt-i]+= np.mean((f1[i:]-f1[i:].mean())*(f2[:-i]-f2[:-i].mean()))
                    corr[dt+i]+= np.mean((f1[:-i]-f1[:-i].mean())*(f2[i:]-f2[i:].mean()))

    return corr


if __name__=="__main__":
    args = parse_args()
    prepared_json = args.json

    pprint("Processing {}".format(prepared_json))
    runner = process(make_config(prepared_json, args))

    runner.align()
    min_freq = 0.003
    weighted_global_average = hasattr(runner, 'tree_leaves')

    # While this could be in a wrapper, it is hopefully more readable this way!
    if runner.config["estimate_mutation_frequencies"]:
        runner.global_frequencies(min_freq, average_global=weighted_global_average)

        if not os.path.exists("processed/rising_mutations/"):
            os.makedirs("processed/rising_mutations/")
        for region in ['AS', 'NA', 'EU']:
           mlist = rising_mutations(runner.mutation_frequencies,
                        runner.mutation_frequency_counts,
                     ['HA1', 'HA2'], region=region, dn=4, offset=2,
                     fname = "processed/rising_mutations/%s_%s_rising_mutations.txt"%("_".join(runner.info["prefix"].split('_')[:-2]), region))


    if runner.config["build_tree"]:
        if weighted_global_average: # subsample alignment
            runner.seqs.aln = MultipleSeqAlignment([v for v in runner.seqs.aln
                                                    if v.name in runner.tree_leaves])
            print("subsampled alignment to %d sequences"%len(runner.seqs.aln))
        else:
            print("using alignment as is, no further subsampling for tree building")

        runner.build_tree()
        runner.timetree_setup_filter_run()
        runner.run_geo_inference()

        # estimate tree frequencies here.
        if runner.config["estimate_tree_frequencies"]:
            pivots = runner.get_pivots_via_spacing()
            runner.estimate_tree_frequencies(pivots=pivots)
            for regionTuple in runner.info["regions"]:
                runner.estimate_tree_frequencies(region=str(regionTuple[0]))

        # titers
        if hasattr(runner, "titers"):
            HI_model(runner)
            if runner.info["lineage"] == "h3n2":
                H3N2_scores(runner, runner.tree.tree, runner.config["epitope_mask"])
            if runner.config["auspice"]["titers_export"]:
                HI_export(runner)

        runner.matchClades(clade_designations[runner.info["lineage"]])

        # Predict fitness.
        if runner.config["annotate_fitness"]:
            fitness_model = runner.annotate_fitness()
            print("Fitness model parameters: %s" % str(zip(fitness_model.predictors, fitness_model.model_params)))
            print("Fitness model deviations: %s" % str(zip(fitness_model.predictors, fitness_model.global_sds)))
            print("Abs clade error: %s" % fitness_model.clade_fit(fitness_model.model_params))
            runner.fitness_model = fitness_model

        if not os.path.exists("processed/recurring_mutations/"):
            os.makedirs("processed/recurring_mutations/")
        recurring_mutations(runner.tree.tree,
                            fname_by_position = "processed/recurring_mutations/%s_recurring_positions.txt"%(runner.info["prefix"]),
                            fname_by_mutation = "processed/recurring_mutations/%s_recurring_mutations.txt"%(runner.info["prefix"]))

        # runner.save_as_nexus()
    runner.auspice_export()
