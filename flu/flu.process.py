from __future__ import print_function
import os, sys, glob
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.process import process
from base.utils import fix_names
from flu_titers import HI_model, HI_export, H3N2_scores
from flu_info import clade_designations
import argparse
import numpy as np
from pprint import pprint
from pdb import set_trace
from Bio.Align import MultipleSeqAlignment

def collect_args():
    parser = argparse.ArgumentParser(description = "Process (prepared) JSON(s)")
    parser.add_argument('-j', '--jsons', default=["all"], nargs='+', type=str, help="prepared JSON(s). \"all\" will do them all. Default = all")
    parser.add_argument('--no_mut_freqs', default=False, action='store_true', help="skip mutation frequencies")
    parser.add_argument('--no_tree_freqs', default=False, action='store_true', help="skip tree (clade) frequencies")
    parser.add_argument('--no_tree', default=False, action='store_true', help="skip tree building")
    parser.add_argument('--clean', default=False, action='store_true', help="clean build (remove previous checkpoints)")
    return parser.parse_args()

def make_config (prepared_json, args):
    return {
        "dir": "flu",
        "in": prepared_json,
        "geo_inference": ['region'],
        "auspice": { ## settings for auspice JSON export
            "extra_attr": ['serum'],
            "color_options": {
                "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
            },
            "controls": {'authors':['authors']},
            "defaults": {'geoResolution': ['region'], 'mapTriplicate': True}
        },
        "titers": {
            "criterium": lambda x: len(x.aa_mutations['HA1']+x.aa_mutations['HA2'])>0,
            "epitope_mask": "metadata/h3n2_epitope_masks.tsv",
            "lam_avi":2.0,
            "lam_pot":0.3,
            "lam_drop":2.0
        },
        "build_tree": not args.no_tree,
        "estimate_mutation_frequencies": not args.no_mut_freqs,
        "estimate_tree_frequencies": not args.no_tree_freqs,
        "clean": args.clean,
        "pivot_spacing": 1.0/12,
        "timetree_options": {
            "Tc": 0.03,
            # "confidence":True,
            # "use_marginal":True
        },
        "newick_tree_options":{
            "raxml": False
        }
    }

def rising_mutations(freqs, genes, region='NA', dn=5, baseline = 0.01, fname='tmp.txt'):
    dx = {}
    for gene in genes:
        for mut,f  in freqs[(region, gene)].iteritems():
            tmp_x = f[-dn:].mean()
            tmp_dx = f[-1] - f[-dn]
            dx[(region, gene, mut[0], mut[1])] = (tmp_x, tmp_dx, tmp_dx/(tmp_x+baseline))

    with open(fname, 'w') as ofile:
        ofile.write("#Frequency change over the last %d month in region %s\n"%(dn, region))
        print("#Frequency change over the last %d month in region %s"%(dn, region))
        for k,v in sorted(dx.items(), key=lambda x:x[1][2])[-10:]:
            ofile.write("%s:%d%s: x=%1.3f, dx=%1.3f, dx/x=%1.3f\n"%(k[1], k[2]+1, k[3], v[0], v[1], v[2]))
            print("%s:%d%s: x=%1.3f, dx=%1.3f, dx/x=%1.3f"%(k[1], k[2]+1, k[3], v[0], v[1], v[2]))
        ofile.write('\n')
    return dx

if __name__=="__main__":
    args = collect_args()
    jsons = glob.glob("prepared/*.json") if "all" in args.jsons else args.jsons

    for prepared_json in jsons:
        pprint("Processing {}".format(prepared_json))
        runner = process(make_config(prepared_json, args))
        runner.align()
        min_freq = 0.01
        # estimate mutation frequencies here.
        # While this could be in a wrapper, it is hopefully more readable this way!
        if runner.config["estimate_mutation_frequencies"]:
            runner.seqs.diversity_statistics()
            include_set = {}
            for prot in ['nuc'] + runner.seqs.translations.keys():
                include_set[prot] = np.where(np.sum(runner.seqs.af[prot][:-2]**2, axis=0)<np.sum(runner.seqs.af[prot][:-2], axis=0)**2-min_freq)[0]

            pivots = runner.get_pivots_via_spacing()
            # runner.estimate_mutation_frequencies(pivots=pivots, min_freq=0.02, inertia=np.exp(-1.0/12), stiffness=0.8*12)
            acronyms = set([x[1] for x in runner.info["regions"] if x[1]!=""])
            region_groups = {str(x):[str(y[0]) for y in runner.info["regions"] if y[1] == x] for x in acronyms}
            pop_sizes = {str(x):np.sum([y[-1] for y in runner.info["regions"] if y[1] == x]) for x in acronyms}
            for region in region_groups.iteritems():
                runner.estimate_mutation_frequencies(pivots=pivots, region=region, min_freq=0.02, include_set=include_set,
                                                     inertia=np.exp(-2.0/12), stiffness=2.0*12)
            total_popsize = np.sum(pop_sizes.values())
            for prot in ['nuc'] + runner.seqs.translations.keys():
                weights = {region: np.array(runner.mutation_frequency_counts[region], dtype = float)
                           for region in acronyms}
                for region in weights:
                    weights[region] = np.maximum(0.1, weights[region]/weights[region].max())
                    weights[region]*=pop_sizes[region]

                total_weight = np.sum([weights[region] for region in acronyms],axis=0)
                gl_freqs, gl_counts, gl_confidence = {}, {}, {}
                all_muts = set()
                for region in acronyms:
                    all_muts.update(runner.mutation_frequencies[(region, prot)].keys())
                for mut in all_muts:
                    gl_freqs[mut] = np.sum([runner.mutation_frequencies[(region, prot)][mut]*weights[region] for region in acronyms
                                            if mut in runner.mutation_frequencies[(region, prot)]], axis=0)/total_weight
                    gl_confidence[mut] = np.sqrt(np.sum([runner.mutation_frequency_confidence[(region, prot)][mut]**2*weights[region]
                                                         for region in acronyms
                                                if mut in runner.mutation_frequencies[(region, prot)]], axis=0)/total_weight)
                gl_counts = np.sum([runner.mutation_frequency_counts[region] for region in acronyms
                                            if mut in runner.mutation_frequencies[(region, prot)]], axis=0)
                runner.mutation_frequencies[("global", prot)] = gl_freqs
                runner.mutation_frequency_counts["global"] = gl_counts
                runner.mutation_frequency_confidence[("global", prot)] = gl_confidence

            for region in ['AS', 'NA']:
                mlist = rising_mutations(runner.mutation_frequencies, ['HA1', 'HA2'], region=region, fname = "rising_mutations/%s_%s_rising_mutations.txt"%(runner.info["prefix"], region))

        if runner.config["build_tree"]:
            if hasattr(runner, 'tree_leaves'): # subsample alignment
                runner.seqs.aln = MultipleSeqAlignment([v for v in runner.seqs.aln
                                                        if v.name in runner.tree_leaves])
                print("subsampled alignment to %d sequences"%len(runner.seqs.aln))
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
            if runner.config["titers"]:
                HI_model(runner)
                H3N2_scores(runner, runner.tree.tree, runner.config["titers"]["epitope_mask"])
                HI_export(runner)

            runner.matchClades(clade_designations[runner.info["lineage"]])

            # runner.save_as_nexus()
        runner.auspice_export()

