from __future__ import print_function
import os, sys
sys.path.append('..') # we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
from base.process import process
import argparse

##### Define clades #####
genotypes = {
'denv1': {'I': [('E', 461, 'V'), ('E', 484, 'L'), ('M', 107, 'T')],
           'II': [('E', 345, 'A'), ('E', 432, 'M'), ('E', 439, 'V')],
           'IV': [('E', 339, 'S'), ('M', 72, 'E'), ('E', 88, 'T')],
           'V': [('E', 297, 'T'), ('NS5', 135, 'M')]},
 'denv2': {'AMERICAN': [('E', 71, 'D'), ('E', 81, 'T'), ('E', 129, 'I')],
           'ASIANAMERICAN': [('E', 491, 'A'), ('M', 15, 'G'), ('M', 39, 'I')],
           'ASIANI': [('E', 484, 'I'), ('NS5', 688, 'I'), ('NS1', 222, 'N')],
           'COSMOPOLITAN': [('E', 71, 'A'), ('E', 149, 'N'), ('E', 462, 'V')],
           'SYLVATIC': [('E', 59, 'F'), ('E', 236, 'M'), ('E', 432, 'V')]},
 'denv3': {'I': [('E', 233, 'K'), ('M', 128, 'F'), ('E', 68, 'V')],
           'II': [('M', 57, 'A'), ('NS5', 750, 'K')],
           'III': [('E', 303, 'T'), ('E', 454, 'V'), ('E', 132, 'Y')],
           'IV': [('E', 22, 'E'), ('E', 50, 'V'), ('E', 62, 'G')]},
 'denv4': {'I': [('E', 494, 'H'), ('NS1', 5, 'A')],
           'II': [('E', 265, 'A'), ('E', 46, 'T'), ('NS1', 246, 'S')],
           'SYLVATIC': [('E', 132, 'V'), ('E', 154, 'S'), ('E', 162, 'T')]},
 'all': {}
 }

for i in ['denv1', 'denv2', 'denv3', 'denv4']:
    for k,v in genotypes[i].iteritems():
        genotypes['all'][i.upper()+'_'+k] = v

regions = ['africa', 'europe', 'north_america', 'china', 'south_asia',
            'japan_korea', 'south_pacific', 'oceania', 'south_america',
            'southeast_asia', 'west_asia']


def titer_model(process, **kwargs):
    '''
    estimate a titer tree and substitution model using titers in titer_fname.
    '''
    from base.titer_model import TreeModel, SubstitutionModel
    ## TREE MODEL
    process.titer_tree = TreeModel(process.tree.tree, process.titers, **kwargs)
    process.titer_tree.prepare(**kwargs) # make training set, find subtree with titer measurements, and make_treegraph
    process.titer_tree.train(**kwargs) # pick longest branch on path between each (test, ref) pair, assign titer drops to this branch
                             # then calculate a cumulative antigenic evolution score for each node
    # add tree attributes to the list of attributes that are saved in intermediate files
    for n in process.tree.tree.find_clades():
        n.attr['cTiter'] = n.cTiter
        n.attr['dTiter'] = n.dTiter

    # SUBSTITUTION MODEL
    process.titer_subs = SubstitutionModel(process.tree.tree, process.titers, **kwargs)
    process.titer_subs.prepare(**kwargs)
    process.titer_subs.train(**kwargs)

    if kwargs['training_fraction'] != 1.0:
        process.titer_tree.validate() #(plot=True, fname='treeModel_%s.png'%lineage)
        process.titer_subs.validate() #(plot=True, fname='subModel_%s.png'%lineage)

    process.config["auspice"]["color_options"]["cTiter"] = {
    "menuItem": "antigenic advance", "type": "continuous", "legendTitle": "log2 titer distance from root", "key": "cTiter", "vmin": "0.0", "vmax": "2.0"
        }


def titer_export(process):
    from base.io_util import write_json
    prefix = process.config["output"]["auspice"]+'/'+process.info["prefix"]+'_'
    if hasattr(process, 'titer_tree'):
        # export the raw titers
        data = process.titer_tree.compile_titers()
        write_json(data, prefix+'titers.json', indent=1)
        # export the tree model (avidities and potencies only)
        tree_model = {'potency':process.titer_tree.compile_potencies(),
                      'avidity':process.titer_tree.compile_virus_effects(),
                      'dTiter':{n.clade:n.dTiter for n in process.tree.tree.find_clades() if n.dTiter>1e-6}}
        write_json(tree_model, prefix+'tree_model.json')
    else:
        print('Tree model not yet trained')

    if hasattr(process, 'titer_tree'):
        # export the substitution model
        titer_subs_model = {'potency':process.titer_subs.compile_potencies(),
                      'avidity':process.titer_subs.compile_virus_effects(),
                      'substitution':process.titer_subs.compile_substitution_effects()}
        write_json(titer_subs_model, prefix+'titer_subs_model.json')
    else:
        print('Substitution model not yet trained')


def collect_args():
    parser = argparse.ArgumentParser(description = "Process (prepared) JSON(s)")
    parser.add_argument('-j', '--jsons', '--json', default=None, nargs='+', type=str, help="Accepts path to prepared JSON(s); overrides -s argument")
    parser.add_argument('-s', '--serotypes', default=["multiple"], nargs='+', type=str, choices=['denv1', 'denv2', 'denv3', 'denv4', 'all', 'multiple'],
    help="Look for prepared JSON(s) like ./prepared/dengue_SEROTYPE.json; 'multiple' will run all five builds. Default='multiple'")
    parser.add_argument('--clean', default=False, action='store_true', help="clean build (remove previous checkpoints)")
    parser.add_argument('--no_mut_freqs', default=False, action='store_true', help="skip mutation frequencies")
    parser.add_argument('--no_tree_freqs', default=False, action='store_true', help="skip tree (clade) frequencies")
    parser.add_argument('--no_titers', default=False, action='store_true', help="skip titer models")
    return parser.parse_args()

def make_config (prepared_json, args):
    return {
        "dir": "dengue",
        "in": prepared_json,
        "geo_inference": ['region'], # what traits to perform this on
        "auspice": { ## settings for auspice JSON export
            "extra_attr": ['serum'],
            "color_options": {
                # "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
                "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
                "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"},
            },
            "controls": {'authors':['authors']},
            "defaults": {'geoResolution': ['region'], 'colorBy': ['region'], 'distanceMeasure': ['div'], 'mapTriplicate': True}
            },

        "timetree_options": {"Tc": False},
        "fit_titer_model": not args.no_titers,
        "titers": {
            "lam_avi":3.0,
            "lam_pot":0.5,
            "lam_drop":1.0,
            "training_fraction":0.9
        },
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

        # estimate mutation frequencies here.
        # if runner.config["estimate_mutation_frequencies"]:
        #     pivots = runner.get_pivots_via_spacing()
        #     runner.estimate_mutation_frequencies(pivots=pivots, min_freq=0.02, inertia=np.exp(-1.0/12), stiffness=0.8*12)

        # estimate tree frequencies here.
        if runner.config["estimate_tree_frequencies"]:
            pivots = runner.get_pivots_via_spacing()
            runner.estimate_tree_frequencies(pivots=pivots)

            for region in ['southeast_asia', 'south_america']: #regions:
                try:
                    runner.estimate_tree_frequencies(region=region)
                except:
                    continue
        # titers
        if runner.config["fit_titer_model"] and runner.config["titers"]:
            titer_model(runner,
                        lam_pot = runner.config['titers']['lam_pot'],
                        lam_avi = runner.config['titers']['lam_avi'],
                        lam_drop = runner.config['titers']['lam_drop'],
                        training_fraction = runner.config['titers']['training_fraction'])
            titer_export(runner)

        runner.matchClades(genotypes[runner.info['lineage']])

        # runner.save_as_nexus()
        runner.auspice_export()
