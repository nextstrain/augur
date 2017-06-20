from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.process import process
import argparse

def collect_args():
    parser = argparse.ArgumentParser(description = "Process (prepared) JSON(s)")
    parser.add_argument('-j', '--jsons', '--json', '-s', '--serotypes', default=["multiple"], nargs='+', type=str, help="Accepts path to prepared JSON(s) or names of serotypes.\n If handed serotype(s) [\"denv1\", \"denv2\", \"denv3\", \"denv4\", or \"all\"], will look for ./prepared/dengue_SEROTYPE.json; \n \"multiple\" will run all five builds. Default = multiple")
    parser.add_argument('--clean', default=False, action='store_true', help="clean build (remove previous checkpoints)")

    ## Not yet working for dengue; turn off by default for now.
    parser.add_argument('--no_mut_freqs', default=True, action='store_true', help="skip mutation frequencies")
    parser.add_argument('--no_tree_freqs', default=True, action='store_true', help="skip tree (clade) frequencies")
    return parser.parse_args()

def make_config (prepared_json, args):
    return {
        "dir": "dengue",
        "in": prepared_json,
        "geo_inference": ['region'], # what traits to perform this on
        "auspice": { ## settings for auspice JSON export
            # "extra_attr": ['serum'],
            "color_options": {
                # "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
                "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
            },
            "controls": {'authors':['authors']},
            "defaults": {'geoResolution': ['region'], 'colorBy': ['region'], 'distanceMeasure': ['div'], 'mapTriplicate': True}
            },

        "timetree_options": {"Tc": False},
        # "titers": {"fname": "../../fauna/data/<LINEAGE>_titers.tsv"},
        "estimate_mutation_frequencies": not args.no_mut_freqs,
        "estimate_tree_frequencies": not args.no_tree_freqs,
        "clean": args.clean,
        "pivot_spacing": 1.0/12, # is this n timepoints per year?
    }

if __name__=="__main__":
    args = collect_args()

    serotypes = ['denv1', 'denv2', 'denv3', 'denv4', 'all']
    if "multiple" in args.jsons: # run all 5 builds from the JSONs in ./prepared/
        args.jsons = serotypes
    for i, j in enumerate(args.jsons): # if passed serotypes instead of paths, validate the corresponding JSONs exist in ./prepared/
        if not os.path.isfile(j):
            assert j in serotypes, "ERROR: %s is not a valid JSON or serotype.\nPass either a valid JSON path name or, to automatically look for ./prepared/dengue_SEROTYPE.json, provide a valid serotype: ['denv1', 'denv2', 'denv3', 'denv4', 'all']"%j

            assert os.path.isfile('./prepared/dengue_%s.json'%j), 'ERROR: no JSON found for serotype %s'%j
            args.jsons[i] = './prepared/dengue_%s.json'%j

    for prepared_json in args.jsons:
        print("Processing %s"%prepared_json)
        runner = process(make_config(prepared_json, args))
        runner.align()
        runner.build_tree()
        runner.timetree_setup_filter_run()
        runner.run_geo_inference()

        ############# WORK IN PROGRESS ############

        # estimate mutation frequencies here.
        # if runner.config["estimate_mutation_frequencies"]:
        #     pivots = runner.get_pivots_via_spacing()
        #     runner.estimate_mutation_frequencies(pivots=pivots, min_freq=0.02, inertia=np.exp(-1.0/12), stiffness=0.8*12)
        #     acronyms = set([x[1] for x in runner.info["regions"] if x[1]!=""])
        #     region_groups = {str(x):[str(y[0]) for y in runner.info["regions"] if y[1] == x] for x in acronyms}
            # for region in region_groups.iteritems():
            #     runner.estimate_mutation_frequencies(region=region, min_freq=0.02, inertia=np.exp(-1.0/12), stiffness=0.8*12)

        # # estimate tree frequencies here.
        # if runner.config["estimate_tree_frequencies"]:
        #     pivots = runner.get_pivots_via_spacing()
        #     runner.estimate_tree_frequencies(pivots=pivots)
        #     for regionTuple in runner.info["regions"]:
        #         runner.estimate_tree_frequencies(region=str(regionTuple[0]))

        # # titers
        # if runner.config["titers"]:
        #     HI_model(runner, )
        #     # H3N2_scores(runner.tree.tree, runner.config["titers"]["epitope_mask"])
        #     HI_export(runner)

        # runner.matchClades(clade_designations[runner.info["lineage"]])

        # runner.save_as_nexus()


        runner.auspice_export()


    # config["in"] = params.json
