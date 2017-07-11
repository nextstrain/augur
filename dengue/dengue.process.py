from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.process import process
import argparse

##### Define clades #####
genotypes = {
                        "denv1":{
                        'I': [('E', 8, 'S'),('E', 37, 'D'),('E', 155, 'S')],
                        'II': [('E', 180, 'T'),('E', 345, 'A')],
                        "V": [('E', 114, 'L'),('E', 161, 'I'),('E', 297, 'T')],
                        "IV": [('E', 37, 'D'),('E', 88, 'T'),]},

                        "denv2": {
                        'ASIAN/AMERICAN': [('E', 203, 'D'), ('E', 491, 'A')],
                        'AMERICAN': [('E', 71, 'D'), ('E', 81, 'T'), ('E', 203, 'D')],
                        'ASIAN-I': [('E', 83, 'K'), ('E', 141, 'V'), ('E', 226, 'K')],
                        'COSMOPOLITAN': [('E', 71, 'A'), ('E', 149, 'N'), ('E', 390, 'S')],
                        'SYLVATIC': [('E', 59, 'F'), ('E', 83, 'V'), ('E', 129, 'V')]},

                        "denv3": {
                        'I': [('E', 68, 'V'),('E', 124, 'S'),('E', 233, 'K')],
                        'II': [('E', 154, 'D')],
                        'III': [('E', 81, 'V'),('E', 132, 'Y'),('E', 171, 'T'),],
                        'IV': [('E', 22, 'E'),('E', 50, 'V'),('E', 62, 'G')]},

                        "denv4": {
                        'I': [('E', 233, 'H'),('E', 329, 'T'),('E', 429, 'L')],
                        'II': [('E', 46, 'T'), ('E', 478, 'T')],
                        'SYLVATIC': [('E', 132, 'V'), ('E', 154, 'S')]},

                        "all": {}
                }

for i in ['denv1', 'denv2', 'denv3', 'denv4']:
    for k,v in genotypes[i].iteritems():
        genotypes['all'][i.upper()+' '+k] = v

def collect_args():
    parser = argparse.ArgumentParser(description = "Process (prepared) JSON(s)")
    parser.add_argument('-j', '--jsons', '--json', default=None, nargs='+', type=str, help="Accepts path to prepared JSON(s); overrides -s argument")
    parser.add_argument('-s', '--serotypes', default=["multiple"], nargs='+', type=str, choices=['denv1', 'denv2', 'denv3', 'denv4', 'all', 'multiple'],
    help="Look for prepared JSON(s) like ./prepared/dengue_SEROTYPE.json; 'multiple' will run all five builds. Default='multiple'")
    parser.add_argument('--clean', default=False, action='store_true', help="clean build (remove previous checkpoints)")
    parser.add_argument('--no_mut_freqs', default=False, action='store_true', help="skip mutation frequencies")
    parser.add_argument('--no_tree_freqs', default=False, action='store_true', help="skip tree (clade) frequencies")
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

    if not args.jsons:
        if 'multiple' in args.serotypes: # "multiple" = run all 5 builds
            args.serotypes = ['denv1', 'denv2', 'denv3', 'denv4', 'all']
        args.jsons = ['./prepared/dengue_%s.json'%s for s in args.serotypes] # Look for ./prepared/dengue_SEROTYPE.json if no file paths given

    for j in args.jsons:            # validate input JSONs exist
        assert os.path.isfile(j)

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

        # estimate tree frequencies here.
        if runner.config["estimate_tree_frequencies"]:
            pivots = runner.get_pivots_via_spacing()
            print('pivots', pivots)
            runner.estimate_tree_frequencies(pivots=pivots)
            for regionTuple in runner.info["regions"]:
                runner.estimate_tree_frequencies(region=str(regionTuple[0]))

        # # titers
        # if runner.config["titers"]:
        #     HI_model(runner, )
        #     # H3N2_scores(runner.tree.tree, runner.config["titers"]["epitope_mask"])
        #     HI_export(runner)

        runner.matchClades(genotypes[runner.info['lineage']])

        runner.save_as_nexus()
        runner.auspice_export()


    # config["in"] = params.json
