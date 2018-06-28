from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import base.process
from base.process import process
import argparse
import numpy as np
from dengue_titers import titer_model, titer_export ## Set up and parameterize the titer model separately for tidiness

##### Define references and metadata #####
sanofi_vaccine_strains = {
    'denv1': 'DENV1/THAILAND/PUO359/1980',
    'denv2': 'DENV2/THAILAND/PUO218/1980',
    'denv3': 'DENV3/THAILAND/PAH88188/1988',
    'denv4': 'DENV4/INDONESIA/S1228/1978',
    'all': None}

regions = ['africa', 'europe', 'north_america', 'china', 'south_asia',
            'japan_korea', 'south_pacific', 'oceania', 'south_america',
            'southeast_asia', 'west_asia']

##### Parse args, set up config #####
def collect_args():
    """Returns a dengue-specific argument parser."""
    parser = base.process.collect_args()

    # parser.add_argument('-j', '--json', default=None, nargs='+', type=str, help="Accepts path to prepared JSON(s); overrides -s argument")
    parser.add_argument('-s', '--serotypes', default=["multiple"], nargs='+', type=str, choices=['denv1', 'denv2', 'denv3', 'denv4', 'all', 'multiple'],
    help="Look for prepared JSON(s) like ./prepared/dengue_SEROTYPE.json; 'multiple' will run all five builds. Default='multiple'")
    parser.add_argument('--no_mut_freqs', default=True, action='store_true', help="skip mutation frequencies")
    parser.add_argument('--no_tree_freqs', default=False, action='store_true', help="skip tree (clade) frequencies")
    parser.add_argument('--no_titers', default=False, action='store_true', help="skip titer models")
    parser.set_defaults(json = None)
    return parser


def make_config (prepared_json, args):
    """
    Configure your analysis here.
    Parsed as a function to enable running multiple builds with one cmd.
    """
    return {
        "dir": "dengue",
        "in": prepared_json,
        "geo_inference": ['region'], # what traits to perform this on; don't run country (too many demes, too few sequences per deme to be reliable)
        "auspice": { ## settings for auspice JSON export
            "extra_attr": ['serum', 'clade', 'dTiter_sanofi'], # keys from tree.tree.clade['attr'] to include in export
            "color_options": { # which traits to color the tree by in auspice; titer colorbys are added in dengue_titers
                "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
                "region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
                "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"},
            },
            "defaults": {'geoResolution': 'region', 'colorBy': 'region', 'distanceMeasure': 'div', 'mapTriplicate': True}
            },

        "timetree_options": {"Tc": False},
        "fit_titer_model": not args.no_titers,
        "titers": { # regularization parameter values and cross-validation fraction
            "lam_avi":0.0,
            "lam_pot":0.5,
            "lam_drop":1.0,
            "training_fraction":0.9,
        },
        "estimate_mutation_frequencies": not args.no_mut_freqs,
        "estimate_tree_frequencies": not args.no_tree_freqs,
        "clean": args.clean,
        "pivot_spacing": 1.0/4, # pivots = time points; 1/N timepoints per year
        "newick_tree_options":{
            "raxml": not args.no_raxml # for dev work only
        }
    }

##### Parse input files/params and run #####
if __name__=="__main__":
    parser = collect_args()
    args = parser.parse_args()

    ### Find the right input files ###
    if args.json: # If provided, a specified JSON path overrides serotype argument
        args.json = [args.json]
    else: # Look for serotype-specific JSONs in the ./prepared/ directory
        if 'multiple' in args.serotypes: # "multiple" = run all 5 builds
            args.serotypes = ['denv1', 'denv2', 'denv3', 'denv4', 'all']
        else:
            args.serotypes = args.serotypes
        args.json = ['./prepared/dengue_%s.json'%s for s in args.serotypes] # Look for ./prepared/dengue_SEROTYPE.json if no file paths given

    for j in args.json:            # validate input JSONs exist
        assert os.path.isfile(j)


    ### Run analyses ###
    for prepared_json in args.json:
        try:
            print("Processing %s"%prepared_json)
            runner = process(make_config(prepared_json, args)) # parse
            runner.align() # run alignment with mafft
            runner.build_tree() # build tree with fasttree -> raxml
            runner.timetree_setup_filter_run() # infer ML ancestral states (geo traits, node dates, mutations)
            runner.run_geo_inference() # run mugration model to infer transmissions

            # estimate mutation frequencies here.
            if runner.config["estimate_mutation_frequencies"]:
                pivots = runner.get_pivots_via_spacing()
                runner.estimate_mutation_frequencies(pivots=pivots, min_freq=0.02, inertia=np.exp(-1.0/12), stiffness=2)

            # estimate tree frequencies here.
            if runner.config["estimate_tree_frequencies"]: # methods @ [ref]
                pivots = runner.get_pivots_via_spacing()
                runner.estimate_tree_frequencies(pivots=pivots, stiffness=2) # stiffness ~= amount of smoothing

                for region in ['southeast_asia', 'south_america']: #regions:
                    try:
                        runner.estimate_tree_frequencies(region=region, stiffness=2)
                    except:
                        continue
            # titers
            if runner.config["fit_titer_model"] and runner.config["titers"]: # methods @ Neher et al., PNAS 2016
                titer_model(runner,
                            lam_pot = runner.config['titers']['lam_pot'],
                            lam_avi = runner.config['titers']['lam_avi'],
                        lam_drop = runner.config['titers']['lam_drop'],
                        training_fraction = runner.config['titers']['training_fraction'],
                        sanofi_strain = sanofi_vaccine_strains[runner.info['lineage']], # vaccine strain for each serotype-specific build
                            plot=False,
                            criterium = lambda node: True) # calculate dTiter for all branches
                        cross_validate=3) # calculate dTiter for all branches
            titer_export(runner)

        ### Export for visualization in auspice
        runner.auspice_export()

        except:
            continue

##### Extra code bank #####
'''
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

# Label named clades based on mutations/genotypes at defining sites
runner.matchClades(genotypes[runner.info['lineage']])
# this is tricky with dengue because the canonical genotypes
# don't really represent the present-day viral diversity.
# I'll get around to redefining these soon-ish hopefully.

### Comparison: force dTiter values to be non-zero only on interserotype brances
    def is_interserotype(node):
        descendents = node.get_terminals()
        serotypes = [k.name.split('/')[0] for k in descendents if 'DENV' in k.name]
        serotypes = [s for s in serotypes if s != 'DENV']
        return len(set(serotypes)) > 1

    interserotype_branches = []
    for node in runner.tree.tree.find_clades():
        if is_interserotype(node):
            interserotype_branches.append(node)
            for child in node.clades:
                interserotype_branches.append(child)
    for node in runner.tree.tree.find_clades():
        if node in interserotype_branches:
            node.interserotype = True
        else:
            node.interserotype = False

    titer_model(runner,
                lam_pot = runner.config['titers']['lam_pot'],
                lam_avi = runner.config['titers']['lam_avi'],
                lam_drop = runner.config['titers']['lam_drop'],
                training_fraction = runner.config['titers']['training_fraction'],
                plot=False,
                criterium = lambda node: node.interserotype == True,
                csv_fname='~/Users/Sidney/Dropbox/dengue/data/titer-model/interserotype-branch-effects/model_predictions.csv')
'''
