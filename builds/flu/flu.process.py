from __future__ import print_function
import os, sys, glob
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import base.process
from base.fitness_model import process_predictor_args
from base.process import process
from base.utils import fix_names
from flu_titers import HI_model, HI_export, H3N2_scores, seasonal_flu_scores
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
            "criterium": lambda x: sum([len(x.aa_mutations[k]) for k in x.aa_mutations])>0,
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

def flatten_json(j):
    nodes = {}
    if "children" in j:
        for c in j["children"]:
            nodes.update(flatten_json(c))
    nodes[j["strain"]] = j
    return nodes


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


def age_distribution(runner):
    import matplotlib
    # important to use a non-interactive backend, otherwise will crash on cluster
    matplotlib.use('PDF')
    import matplotlib.pyplot as plt

    fs=16
    bins = np.arange(0,100,10)
    bc = 0.5*(bins[1:]+bins[:-1])
    plt.figure()
    plt.title(runner.info["prefix"])
    with open('processed/%s_age_distributions.txt'%runner.info["prefix"], 'w') as ofile:
        ofile.write("\t".join(map(str, ['region', 'total'] + ["%d-%d"%(bins[i], bins[i+1]) for i in range(len(bins)-1)]))+'\n')
        for region in ['total', 'africa','china','europe','japan_korea','north_america','oceania','south_america','south_asia','southeast_asia','west_asia']:
            y,x = np.histogram([n.attr['age'] for n in runner.tree.tree.get_terminals()
                                if n.attr['age']!='unknown' and (n.attr['region']==region or region=='total')], bins=bins)
            total = np.sum(y)
            y = np.array(y, dtype=float)/total
            ofile.write("\t".join(map(str, [region, total]+list(y)))+'\n')
            plt.plot(bc, y, label=region, lw=3 if region=='total' else 2)

    plt.legend(fontsize=fs*0.8, ncol=2)
    plt.ylabel('age distribution', fontsize=fs)
    plt.xlabel('age', fontsize=fs)
    plt.ylim([0,0.7])
    plt.savefig('processed/%s_age_distributions.png'%runner.info["prefix"])
    plt.close()


def plot_titers(titer_model, titers, fname=None, title=None):
    from collections import defaultdict
    import matplotlib
    # important to use a non-interactive backend, otherwise will crash on cluster
    matplotlib.use('PDF')
    import matplotlib.pyplot as plt

    symb = ['o', 's', 'd', 'v', '<', '>', '+']
    cols = ['C'+str(i) for i in range(10)]
    fs =16
    grouped_titers = defaultdict(list)
    for (test, (ref, serum)), val in titers.items():
        if test not in titer_model.node_lookup:
            continue
        date = titer_model.node_lookup[test].attr["num_date"]
        date = int(date*2)/2.
        if ref!=test:
            if np.isscalar(val):
                grouped_titers[(ref, serum, date)].append(val)
            else:
                grouped_titers[(ref, serum, date)].append(np.exp(np.mean(np.log(val))))


    titer_means = defaultdict(list)
    for (ref, serum, date), val in grouped_titers.items():
        d = [date, np.mean(val), np.std(val), len(val)]
        titer_means[(ref, serum)].append(d)

    for k in titer_means:
        titer_means[k] = np.array(sorted(titer_means[k], key=lambda x:x[0]))

    plt.figure(figsize=(16,12))
    if title is not None:
        plt.title(title)
    ii = 0
    curves =  [k for k, val in titer_means.items() if np.sum(val[-3:,-1])>100]
    n_curves = len(curves)
    virus = list(set([k[0] for k in curves]))

    for k in sorted(titer_means.keys()):
        val = titer_means[k]
        if np.sum(val[-3:,-1])>100:
            sub_val = val[val[:,-1]>10]
            c= cols[virus.index(k[0])%len(cols)]
            plt.errorbar(sub_val[:,0]+(ii-0.5*n_curves)*0.01, sub_val[:,1], sub_val[:,2],
                        label=k[0]+', '+k[1], lw=2, c=c,
                        marker=symb[ii%len(symb)], markersize=10)
            ii+=1
    plt.legend(ncol=2, fontsize=fs*0.7)
    plt.ylabel('titer', fontsize=fs)
    plt.xlabel('year', fontsize=fs)

    if fname is not None:
        plt.savefig(fname)


if __name__=="__main__":
    args = parse_args()
    prepared_json = args.json

    print("Processing {}".format(prepared_json))
    runner = process(make_config(prepared_json, args))

    # this should be in the json...
    segment = "ha"
    if "_na_" in prepared_json:
        segment = "na"
    runner.segment = segment

    runner.align()
    min_freq = 0.003
    weighted_global_average = hasattr(runner, 'tree_leaves')

    # While this could be in a wrapper, it is hopefully more readable this way!
    if runner.config["estimate_mutation_frequencies"]:
        runner.global_frequencies(min_freq, average_global=weighted_global_average)

        # so far do this only for HA
        if segment=='ha':
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
        seasonal_flu_scores(runner, runner.tree.tree)
        if hasattr(runner, "titers"):
            HI_model(runner)
            if runner.info["lineage"] == "h3n2":
                H3N2_scores(runner, runner.tree.tree, runner.config["epitope_mask"])
            if runner.config["auspice"]["titers_export"]:
                HI_export(runner)
                plot_titers(runner.HI_subs, runner.HI_subs.titers.titers,
                            fname='processed/%s_raw_titers.png'%runner.info["prefix"],
                            title = runner.info["prefix"])
                plot_titers(runner.HI_subs, runner.HI_subs.titers.titers_normalized,
                            fname='processed/%s_normalized_titers.png'%runner.info["prefix"],
                            title = runner.info["prefix"])

        # outputs figures and tables of age distributions
        age_distribution(runner)
        # ignore fitness for NA.
        if segment=='ha':
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
    if segment=="na":
        import json
        ha_tree_json_fname = os.path.join(runner.config["output"]["auspice"], runner.info["prefix"]) + "_tree.json"
        ha_tree_json_fname = ha_tree_json_fname.replace("_na", "_ha")
        if os.path.isfile(ha_tree_json_fname):      # confirm file exists
            with open(ha_tree_json_fname) as jfile:
                ha_tree_json = json.load(jfile)
            ha_tree_flat = flatten_json(ha_tree_json)

            for n in runner.tree.tree.get_terminals():
                if n.name in ha_tree_flat:
                    if "named_clades" in ha_tree_flat[n.name]["attr"]:
                        n.attr["named_clades"] = ha_tree_flat[n.name]["attr"]["named_clades"]
                else:
                    n.attr["named_clades"] = ["unassigned"]

    runner.auspice_export()
