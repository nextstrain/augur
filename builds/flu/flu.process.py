from __future__ import print_function
import os, sys, glob
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import base.process
from base.fitness_model import process_predictor_args
from base.process import process
from base.utils import fix_names
from base.io_util import write_json
from flu_titers import HI_model, HI_export, vaccine_distance
from scores import calculate_sequence_scores, calculate_metadata_scores, calculate_phylogenetic_scores
from flu_info import clade_designations, lineage_to_epitope_mask, lineage_to_glyc_mask, resolution_to_pivot_spacing, vaccine_choices
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
    parser.add_argument('--pivot_spacing', type=float, help="month per pivot")
    parser.add_argument('--titers_export', default=False, action='store_true', help="export titers.json file")
    parser.add_argument('--annotate_fitness', default=False, action='store_true', help="run fitness prediction model and annotate fitnesses on tree nodes")
    parser.add_argument('--predictors', default=['cTiter'], nargs='+', help="attributes to use as fitness model predictors")
    parser.add_argument('--predictors_params', type=float, nargs='+', help="precalculated fitness model parameters for each of the given predictors")
    parser.add_argument('--predictors_sds', type=float, nargs='+', help="precalculated global standard deviations for each of the given predictors")
    parser.add_argument('--epitope_mask_version', help="name of the epitope mask that defines epitope mutations")
    parser.add_argument('--tolerance_mask_version', help="name of the tolerance mask that defines non-epitope mutations")
    parser.add_argument('--glyc_mask_version', help="name of the mask that defines putative glycosylation sites")

    parser.set_defaults(
        json="prepared/flu.json"
    )

    return parser.parse_args()

def make_config(prepared_json, args):
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
            "panels": ['tree', 'entropy', 'frequencies'],
            "extra_attr": ['serum'],
            "color_options": {
                "region":{"menuItem":"region", "legendTitle":"Region", "key":"region", "type":"discrete"},
                "clade_membership": {"menuItem": "clade", "legendTitle": "Clade", "key": "clade_membership", "type": "discrete"},
            },
            "controls": {'authors':['authors']},
            "defaults": {'colorBy': 'clade_membership',
                'geoResolution': 'region',
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
        "ha_masks": "metadata/ha_masks.tsv",
        "epitope_mask_version": args.epitope_mask_version,
        "tolerance_mask_version": args.tolerance_mask_version,
        "glyc_mask_version": args.glyc_mask_version,
        "annotate_fitness": args.annotate_fitness,
        "predictors": predictors,
        "clean": args.clean,
        "pivot_spacing": args.pivot_spacing,
        "timetree_options": {
            "Tc": 0.03,
            # "confidence":True,
            # "use_marginal":True
        },
        "newick_tree_options":{
            "method": args.tree_method
        }
    }

# set defaults when command line parameter is None based on lineage, segment and resolution
def set_config_defaults(runner):
    if runner.config["epitope_mask_version"] is None and runner.info["lineage"] in lineage_to_epitope_mask:
        runner.config["epitope_mask_version"] = lineage_to_epitope_mask[runner.info["lineage"]]
    if runner.config["glyc_mask_version"] is None and runner.info["lineage"] in lineage_to_glyc_mask:
        runner.config["glyc_mask_version"] = lineage_to_glyc_mask[runner.info["lineage"]]
    if runner.config["pivot_spacing"] is None and runner.info["resolution"] in resolution_to_pivot_spacing:
        runner.config["pivot_spacing"] = resolution_to_pivot_spacing[runner.info["resolution"]]

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
    matplotlib.use('agg')
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


def mean_func(x, mean='arithmetric'):
    if mean=='arithmetric':
        return np.mean(x)
    elif mean=='geometric':
        return np.exp(np.mean(np.log(x)))

def plot_titers(titer_model, titers, fname=None, title=None, mean='geometric'):
    from collections import defaultdict
    import matplotlib
    # important to use a non-interactive backend, otherwise will crash on cluster
    matplotlib.use('agg')
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
                grouped_titers[(ref, serum, date)].append(mean_func(val, mean=mean))


    titer_means = defaultdict(list)
    for (ref, serum, date), val in grouped_titers.items():
        d = [date, mean_func(val, mean=mean), np.std(val), len(val)]
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


def plot_titer_matrix(titer_model, titers, clades=None, fname=None, title=None, mean='arithmetric', normalized=False, potency=False):
    from collections import defaultdict
    import matplotlib
    # important to use a non-interactive backend, otherwise will crash on cluster
    matplotlib.use('agg')
    import matplotlib.pyplot as plt

    symb = ['o', 's', 'd', 'v', '<', '>', '+']
    cols = ['C'+str(i) for i in range(10)]
    fs = 16
    grouped_titers = defaultdict(list)
    autologous = defaultdict(list)
    for (test, (ref, serum)), val in titers.items():
        if test not in titer_model.node_lookup:
            continue
        node = titer_model.node_lookup[test]
        date = node.attr["num_date"]
        date = int(date*2)/2.
        if "named_clades" in node.attr:
            clade = '_'.join(node.attr["named_clades"])
        else:
            clade = 'unassigned'
        if ref!=test:
            if date>=2016:
                if np.isscalar(val):
                    grouped_titers[(ref, serum, clade)].append(val)
                else:
                    grouped_titers[(ref, serum, clade)].append(mean_func(val, mean=mean))
        else:
            autologous[(ref, serum)].append(mean_func(val, mean=mean))

    titer_means = defaultdict(list)
    ntiters = defaultdict(int)
    for k, val in grouped_titers.items():
        d = [date, mean_func(val, mean=mean), np.std(val), len(val)]
        ntiters[(k[0], k[1])]+=d[-1]
        titer_means[k] = d

    for k, val in autologous.items():
        autologous[k] = mean_func(val, mean=mean)

    titer_matrix = []
    sera_with_counts = sorted(ntiters.items(), key=lambda x:x[1], reverse=True)
    rows = []

    sorted_sera_with_counts = []
    # sort sera according to clades
    for clade in clades:
        for serum,count in sera_with_counts:
            serum_strain = serum[0]
            serum_clade = 'unassigned'
            if serum_strain in titer_model.node_lookup:
                node = titer_model.node_lookup[serum_strain]
                if "named_clades" in node.attr:
                    serum_clade = '_'.join(node.attr["named_clades"])
            if clade == serum_clade:
                sorted_sera_with_counts.append((serum,count))

    for serum,count in sorted_sera_with_counts:
        if count>50:
            serum_clade = 'unassigned'
            if serum[0] in titer_model.node_lookup:
                node = titer_model.node_lookup[serum[0]]
                if "named_clades" in node.attr:
                    serum_clade = '_'.join(node.attr["named_clades"])
            tmp_autologous = 640
            if autologous[serum]:
                tmp_autologous = autologous[serum]
            rows.append('\n'.join(serum)+'\n'+serum_clade)
            tmp = []
            for clade in clades:
                k = (serum[0], serum[1], clade)
                meanvalue = np.nan
                if k in titer_means and titer_means[k][-1]>5:
                    if normalized:
                        meanvalue = titer_means[(k)][1]
                    else:
                        if mean=='geometric':
                            meanvalue = -np.log2(titer_means[(k)][1]/tmp_autologous)
                        else:
                            meanvalue = titer_means[(k)][1]-tmp_autologous
                    if potency:
                        meanvalue -= titer_model.serum_potency[serum]
                tmp.append(meanvalue)
            titer_matrix.append(tmp)

    titer_matrix = np.array(titer_matrix)

    if len(rows) > 0:
        import seaborn as sns
        plt.figure(figsize=(12,9))
        if title is not None:
            plt.title(title)
        cmap = sns.cubehelix_palette(start=2.6, rot=.1, as_cmap=True)
        sns.heatmap(titer_matrix, xticklabels=clades, yticklabels=rows,
                    annot=True, cmap=cmap, vmin=0, vmax=4)
        plt.yticks(rotation=0)
        plt.xticks(rotation=30)
        plt.tight_layout()

        if fname is not None:
            plt.savefig(fname)
            plt.close()

# only use with normalized titers (so arithmetric mean)
def plot_titer_matrix_grouped(titer_model, titers, virus_clades=None, serum_clades=None, fname=None, title=None, potency=False, minDate=2016, minMeasurements=30, minRows=2, maxSeraPerClade=5):
    from collections import defaultdict
    import matplotlib
    # important to use a non-interactive backend, otherwise will crash on cluster
    matplotlib.use('agg')
    import matplotlib.pyplot as plt

    symb = ['o', 's', 'd', 'v', '<', '>', '+']
    cols = ['C'+str(i) for i in range(10)]
    fs = 16
    grouped_titers = defaultdict(list)

    autologous_exists = {}
    for (test, (ref, serum)), val in titers.items():
        autologous_exists[ref] = False

    for (test, (ref, serum)), val in titers.items():
        if potency:
            val -= titer_model.serum_potency[(ref, serum)]
        if test not in titer_model.node_lookup:
            continue
        node = titer_model.node_lookup[test]
        date = node.attr["num_date"]
        date = int(date*2)/2.
        if "named_clades" in node.attr:
            clade = '_'.join(node.attr["named_clades"])
        else:
            clade = 'unassigned'
        if ref!=test:
            if date >= minDate and "-egg" not in test: # only keep measurements to cell antigens, egg sera are already normalized
                if np.isscalar(val):
                    grouped_titers[(ref, clade)].append(val)
                else:
                    grouped_titers[(ref, clade)].append(np.mean(val))
        else:
            autologous_exists[ref] = True

    titer_means = defaultdict(list)
    ntiters = defaultdict(int)
    for (serum_strain, clade), val in grouped_titers.items():
        ntiters[serum_strain] += len(val)
        titer_means[(serum_strain, clade)] = [np.mean(val), len(val)]
    sorted_ntiters = sorted(ntiters.items(), key=lambda x:x[1], reverse=True)

    titer_matrix = []
    rows = []
    sera_with_counts = []

    # sort sera according to clade, then by count, this only sorts, keeps all data
    for clade in serum_clades:
        for serum,count in sorted_ntiters:
            serum_strain = serum
            serum_clade = 'unassigned'
            if serum_strain in titer_model.node_lookup:
                node = titer_model.node_lookup[serum_strain]
                if "named_clades" in node.attr:
                    serum_clade = '_'.join(node.attr["named_clades"])
            if clade == serum_clade:
                sera_with_counts.append((serum,count))

    # walk through and keep good bins, taking their average
    count_for_clade = {}
    for clade in serum_clades:
        count_for_clade[clade] = 0

    for serum,count in sera_with_counts:
        if count > minMeasurements and autologous_exists[serum]:
            serum_clade = 'unassigned'
            if serum in titer_model.node_lookup:
                node = titer_model.node_lookup[serum]
                if "named_clades" in node.attr:
                    serum_clade = '_'.join(node.attr["named_clades"])
            tmp = []
            good_rows = 0
            for clade in virus_clades:
                k = (serum, clade)
                meanvalue = np.nan
                if k in titer_means and titer_means[k][-1]>5:
                    meanvalue = titer_means[k][0]
                    good_rows += 1
                tmp.append(meanvalue)
            if good_rows >= minRows and count_for_clade[serum_clade] < maxSeraPerClade:
                rows.append(serum+'\n'+serum_clade)
                titer_matrix.append(tmp)
                count_for_clade[serum_clade] += 1

    titer_matrix = np.array(titer_matrix)

    if len(rows) > 0:
        import seaborn as sns
        plt.figure(figsize=(7, 0.6*len(rows)+1))
        if title is not None:
            title = title.replace("flu_", "").replace("_ha_", "_").replace("_2y_", "_")
            title = title.replace("h3n2", "H3N2").replace("h1n1pdm", "H1N1pdm").replace("vic", "Vic").replace("yam", "Yam")
            title = title.replace("who", "WHO").replace("cdc", "CDC").replace("crick", "Crick").replace("niid", "NIID").replace("vidrl", "VIDRL")
            title = title.replace("hi", "HI").replace("fra", "FRA")
            title = title.replace("_", " ")
            plt.title(title)
        cmap = sns.cubehelix_palette(start=2.6, rot=.1, as_cmap=True)
        ax = sns.heatmap(titer_matrix, xticklabels=virus_clades, yticklabels=rows,
                    annot=True, fmt='2.1f', cmap=cmap, vmin=0, vmax=4)
        plt.yticks(rotation=0)
        plt.xticks(rotation=30)
        cbar = ax.collections[0].colorbar
        cbar.set_ticks([0, 1, 2, 3, 4])
        plt.tight_layout()

        if fname is not None:
            plt.savefig(fname, dpi=200)
            plt.close()


if __name__=="__main__":
    args = parse_args()
    prepared_json = args.json

    print("Processing {}".format(prepared_json))
    runner = process(make_config(prepared_json, args))

    # set defaults
    set_config_defaults(runner)

    runner.align()
    min_freq = 0.003
    weighted_global_average = hasattr(runner, 'tree_leaves')

    # While this could be in a wrapper, it is hopefully more readable this way!
    if runner.config["estimate_mutation_frequencies"]:
        runner.global_frequencies(min_freq, average_global=weighted_global_average)

        # so far do this only for HA
        if runner.info["segment"] in ['ha', 'na']:
            genes_by_segment = {'ha':['HA1', 'HA2'], 'na':['NA']}
            if not os.path.exists("processed/rising_mutations/"):
                os.makedirs("processed/rising_mutations/")
            for region in ['AS', 'NA', 'EU']:
               mlist = rising_mutations(runner.mutation_frequencies,
                            runner.mutation_frequency_counts,
                         genes_by_segment[runner.info["segment"]], region=region, dn=4, offset=2,
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

        if runner.info["segment"]=='ha':
            if runner.info["lineage"]=='h3n2':
                clades = ['3c2.A', 'A1', 'A1b/135K', 'A2', 'A3']
                virus_clades = ['A1', 'A1a', 'A1b/135K', 'A1b/135N', 'A2', 'A3']
                serum_clades = ['3c2.A', 'A1', 'A1a', 'A1b', 'A1b/135K', 'A1b/135N', 'A2', 'A3']
            elif runner.info["lineage"]=='h1n1pdm':
                clades = ['6b.1', '6b.2', '164T']
                virus_clades = clades
                serum_clades = clades
            elif runner.info["lineage"]=='vic':
                clades = ['1A', '117V', 'DV']
                virus_clades = clades
                serum_clades = clades
            elif runner.info["lineage"]=='yam':
                clades = ['2', '3', '172Q']
                virus_clades = clades
                serum_clades = clades
            else:
                clades = clade_designations[runner.info["lineage"]].keys()
                virus_clades = clades
                serum_clades = clades
            runner.matchClades(clade_designations[runner.info["lineage"]])

        if runner.info["segment"] in ['ha', 'na']:
            if not os.path.exists("processed/recurring_mutations/"):
                os.makedirs("processed/recurring_mutations/")
            recurring_mutations(runner.tree.tree,
                                fname_by_position = "processed/recurring_mutations/%s_recurring_positions.txt"%(runner.info["prefix"]),
                                fname_by_mutation = "processed/recurring_mutations/%s_recurring_mutations.txt"%(runner.info["prefix"]))

        # runner.save_as_nexus()
        calculate_metadata_scores(runner.tree.tree)
        assert "age" in runner.tree.tree.root.attr, "age not annotated"

        # outputs figures and tables of age distributions
        age_distribution(runner)

        # runner.config["auspice"]["color_options"]["age"] = {
        #     "menuItem": "average host age in clade",
        #     "type": "continuous",
        #     "legendTitle": "Avg host age in clade",
        #     "key": "age"
        # }
        # runner.config["auspice"]["color_options"]["gender"] = {
        #     "menuItem": "average host gender in clade",
        #     "type": "continuous",
        #     "legendTitle": "Avg host gender in clade",
        #     "key": "num_gender"
        # }

        if "LBI_params" in runner.info:
            calculate_phylogenetic_scores(
                runner.tree.tree,
                tau=runner.info["LBI_params"]["tau"],
                time_window=runner.info["LBI_params"]["time_window"]
            )
            assert "lbi" in runner.tree.tree.root.attr, "LBI not annotated"

            runner.config["auspice"]["color_options"]["lbi"] = {
                "menuItem": "local branching index",
                "type": "continuous",
                "legendTitle": "local branching index",
                "key": "lbi",
                "vmin": 0,
                "vmax": 0.7
            }

        if runner.info["segment"]=='ha' and runner.info["lineage"] in ["h3n2", "h1n1pdm"]:

            print("Calculating scores with epitope mask '%s' and glycosylation mask '%s'." % (runner.config["epitope_mask_version"], runner.config["glyc_mask_version"]))
            calculate_sequence_scores(
                runner.tree.tree,
                runner.config["ha_masks"],
                runner.info["lineage"],
                runner.info["segment"],
                epitope_mask_version=runner.config["epitope_mask_version"],
                glyc_mask_version=runner.config["glyc_mask_version"]
            )
            assert "ep" in runner.tree.tree.root.attr, "epitope mutations not annotated"
            assert "ne" in runner.tree.tree.root.attr, "non-epitope mutations not annotated"
            assert "glyc" in runner.tree.tree.root.attr, "glycosylation not annotated"

            if runner.info["lineage"] == "h3n2":
                assert "rb" in runner.tree.tree.root.attr, "rbs mutations not annotated"

            # Define color options for sequence score annotations.
            runner.config["auspice"]["color_options"]["ep"] = {
                "menuItem": "epitope mutations",
                "type": "continuous",
                "legendTitle": "Epitope mutations",
                "key": "ep"
            }
            runner.config["auspice"]["color_options"]["ne"] = {
                "menuItem": "non-epitope mutations",
                "type": "continuous",
                "legendTitle": "Non-epitope mutations",
                "key": "ne"
            }

            # runner.config["auspice"]["color_options"]["glyc"] = {
            #     "menuItem": "potential glycosylation sites",
            #     "type": "continuous",
            # "legendTitle": "Pot. glycosylation count",
            #     "key": "glyc"
            # }

            if runner.info["lineage"]=='h3n2':
                runner.config["auspice"]["color_options"]["rb"] = {
                    "menuItem": "receptor binding mutations",
                    "type": "continuous",
                    "legendTitle": "Receptor binding mutations",
                    "key": "rb"
                }

        # titers
        if hasattr(runner, "titers") and runner.info["segment"] == "ha":
            HI_model(runner)

            if runner.config["auspice"]["titers_export"]:
                HI_export(runner)
                vaccine_distance_json = vaccine_distance(titer_tree = runner.tree.tree,
                                                         vaccine_strains = vaccine_choices[runner.info['lineage']],
                                                         attributes=['dTiter', 'dTiterSub'])
                write_json(vaccine_distance_json, os.path.join(runner.config["output"]["auspice"], runner.info["prefix"])+'_vaccine_dist.json')

                if runner.info["segment"]=='ha':
                    plot_titers(runner.HI_subs, runner.HI_subs.titers.titers,
                                fname='processed/%s_raw_titers.png'%runner.info["prefix"],
                                title = runner.info["prefix"], mean='geometric')
                    plot_titers(runner.HI_subs, runner.HI_subs.titers.titers_normalized,
                                fname='processed/%s_normalized_titers.png'%runner.info["prefix"],
                                title = runner.info["prefix"], mean='arithmetric')
                    plot_titer_matrix(runner.HI_subs, runner.HI_subs.titers.titers,
                                fname='processed/%s_raw_titer_matrix.png'%runner.info["prefix"],
                                title = runner.info["prefix"], mean='geometric', clades=clades, normalized=False, potency=False)
                    plot_titer_matrix(runner.HI_subs, runner.HI_subs.titers.titers_normalized,
                                fname='processed/%s_normalized_titer_matrix.png'%runner.info["prefix"],
                                title = runner.info["prefix"], mean='arithmetric', clades=clades, normalized=True, potency=False)
                    plot_titer_matrix(runner.HI_subs, runner.HI_subs.titers.titers_normalized,
                                fname='processed/%s_normalized_with_potency_titer_matrix.png'%runner.info["prefix"],
                                title = runner.info["prefix"], mean='arithmetric', clades=clades, normalized=True, potency=True)
                    plot_titer_matrix_grouped(runner.HI_subs, runner.HI_subs.titers.titers_normalized,
                                fname='processed/%s_grouped_titer_matrix.png'%runner.info["prefix"],
                                title = runner.info["prefix"], virus_clades=virus_clades, serum_clades=serum_clades, potency=False)
                    plot_titer_matrix_grouped(runner.HI_subs, runner.HI_subs.titers.titers_normalized,
                                fname='processed/%s_grouped_with_potency_titer_matrix.png'%runner.info["prefix"],
                                title = runner.info["prefix"], virus_clades=virus_clades, serum_clades=serum_clades, potency=True)

    if runner.info["segment"] == "na":
        import json
        ha_tree_json_fname = os.path.join(runner.config["output"]["auspice"], runner.info["prefix"]) + "_tree.json"
        ha_tree_json_fname = ha_tree_json_fname.replace("_na", "_ha")
        if os.path.isfile(ha_tree_json_fname):      # confirm file exists
            with open(ha_tree_json_fname) as jfile:
                ha_tree_json = json.load(jfile)
            ha_tree_flat = flatten_json(ha_tree_json)
            for node in runner.tree.tree.find_clades():
                node.attr['clade_membership'] = 'unassigned'
            for n in runner.tree.tree.get_terminals():
                if n.name in ha_tree_flat:
                    if "clade_membership" in ha_tree_flat[n.name]["attr"]:
                        n.attr["clade_membership"] = ha_tree_flat[n.name]["attr"]["clade_membership"]

    # Predict fitness for HA after all other scores and annotations have completed.
    if runner.info["segment"] == 'ha' and runner.config["annotate_fitness"]:
        fitness_model = runner.annotate_fitness()
        if fitness_model is not None:
            print("Fitness model parameters: %s" % str(zip(fitness_model.predictors, fitness_model.model_params)))
            print("Fitness model deviations: %s" % str(zip(fitness_model.predictors, fitness_model.global_sds)))
            print("Abs clade error: %s" % fitness_model.clade_fit(fitness_model.model_params))
            runner.fitness_model = fitness_model

            runner.config["auspice"]["color_options"]["fitness"] = {
                "menuItem": "fitness",
                "type": "continuous",
                "legendTitle": "Fitness",
                "key": "fitness"
            }

            runner.config["auspice"]["color_options"]["predicted_freq"] = {
                "menuItem": "predicted_freq",
                "type": "continuous",
                "legendTitle": "Predicted frequency",
                "key": "predicted_freq"
            }

    runner.auspice_export()
