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

def collect_args():
    parser = argparse.ArgumentParser(description = "Process (prepared) JSON(s)")
    parser.add_argument('-j', '--jsons', default=["all"], nargs='+', type=str, help="prepared JSON(s). \"all\" will do them all. Default = all")
    parser.add_argument('--no_mut_freqs', default=False, action='store_true', help="skip mutation frequencies")
    parser.add_argument('--no_tree_freqs', default=False, action='store_true', help="skip tree (clade) frequencies")
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
            "fname": "../../fauna/data/<LINEAGE>_crick_hi_cell_titers.tsv",
            "criterium": lambda x: len(x.aa_mutations['HA1']+x.aa_mutations['HA2'])>0,
            "epitope_mask": "metadata/h3n2_epitope_masks.tsv",
            "lam_avi":2.0,
            "lam_pot":0.3,
            "lam_drop":2.0
        },
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
            "raxml": True
        }
    }

if __name__=="__main__":
    args = collect_args()
    jsons = glob.glob("prepared/*.json") if "all" in args.jsons else args.jsons

    for prepared_json in jsons:
        pprint("Processing {}".format(prepared_json))
        runner = process(make_config(prepared_json, args))
        runner.align()

        # estimate mutation frequencies here.
        # While this could be in a wrapper, it is hopefully more readable this way!
        if runner.config["estimate_mutation_frequencies"]:
            pivots = runner.get_pivots_via_spacing()
            runner.estimate_mutation_frequencies(pivots=pivots, min_freq=0.02, inertia=np.exp(-1.0/12), stiffness=0.8*12)
            acronyms = set([x[1] for x in runner.info["regions"] if x[1]!=""])
            region_groups = {str(x):[str(y[0]) for y in runner.info["regions"] if y[1] == x] for x in acronyms}
            for region in region_groups.iteritems():
                runner.estimate_mutation_frequencies(region=region, min_freq=0.02, inertia=np.exp(-1.0/12), stiffness=0.8*12)

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
