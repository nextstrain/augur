import os
from base.io_util import write_json

def HI_model(process):
    '''
    estimate a tree and substitution model using titers.
    '''
    from base.titer_model import TreeModel, SubstitutionModel

    ## define the kwargs explicitly
    kwargs = process.config["titers"]
    # if "criterium" in process.config["titers"]:
    #     kwargs["criterium"] = process.config["titers"]["criterium"]

    ## TREE MODEL
    process.HI_tree = TreeModel(process.tree.tree, process.titers, **kwargs)
    process.HI_tree.prepare(**kwargs)
    process.HI_tree.train(**kwargs)
    # add tree attributes to the list of attributes that are saved in intermediate files
    for n in process.tree.tree.find_clades():
        n.attr['cTiter'] = n.cTiter
        n.attr['dTiter'] = n.dTiter
        # print("cumulative: {} delta: {}".format(n.cTiter, n.dTiter))
    process.config["auspice"]["color_options"]["cTiterSub"] = {
        "menuItem": "antigenic advance", "type": "continuous", "legendTitle": "Antigenic advance", "key": "cTiterSub"
    }

    # SUBSTITUTION MODEL
    process.HI_subs = SubstitutionModel(process.tree.tree, process.titers, **kwargs)
    process.HI_subs.prepare(**kwargs)
    process.HI_subs.train(**kwargs)

    for node in process.tree.tree.find_clades():
        dTiterSub = 0

        if hasattr(node, "aa_muts"):
            for gene, mutations in node.aa_muts.iteritems():
                for mutation in mutations:
                    dTiterSub += process.HI_subs.substitution_effect.get((gene, mutation), 0)

        node.dTiterSub = dTiterSub

        if node.up is not None:
            node.cTiterSub = node.up.cTiterSub + dTiterSub
        else:
            node.cTiterSub = 0

        node.attr["cTiterSub"] = node.cTiterSub
        node.attr["dTiterSub"] = node.dTiterSub

def sum_along_path(path, attribute):
    '''For each tip/node along a path in the tree, pull the value of `attribute`. Return the sum of the attribute values.'''
    if len(path) == 1: # start == end
        return 0.
    try:
        return sum([k.attr[attribute] for k in path]) # assume it's in the .attr dictionary
    except KeyError:
        return sum([getattr(k, attribute) for k in path]) # if that fails, try accessing it directly

def vaccine_distance(titer_tree, vaccine_strains, attributes=['dTiter', 'dTiterSub']):
    '''
    For each tip in the tree, trace the path between the tip and each vaccine strain.
    Sum up the requested attributes along this path.
    Also pull the tip's strain name, numdate, and region.

    Returns a list of dictionaries, which is compatible with write_json
    e.g.,
        [
      {
        "strain": "A/Incheon/677/2006",
        "num_date": 2006.867,
        "region": "japan_korea",
        "ep": {
          "A/Perth/16/2009": 4,
          "A/Victoria/361/2011": 6
        },
        "cTiter": {
          "A/Perth/16/2009": 2.3,
          "A/Victoria/361/2011": 4.2
        },
        "cTiterSub": {
          "A/Perth/16/2009": 1.9,
          "A/Victoria/361/2011": 3.7
        }
      }
    ]
    '''

    tips = [ k for k in titer_tree.get_terminals() ] # all tips in the tree
    vaccine_tips = { k.name: k for k in tips if k.name in vaccine_strains } # { strain_name: biophylo object}

    distance_to_vaccine = []

    for k in tips: # pull basic strain data
        tip_data = { 'strain': k.name,
                     'num_date': k.attr['num_date'],
                     'region': k.attr['region'] }

        # compute the paths between this tip and each vaccine strain once
        paths_to_vaccines = { vaccine_strain: titer_tree.trace(k, vaccine_tip) for vaccine_strain, vaccine_tip in vaccine_tips.items() }

        for a in attributes: # for each attribute of interest, pull the precomputed path and sum up the attribute along it
            summed_attribute = { vaccine_strain: sum_along_path(paths_to_vaccines[vaccine_strain], a) for vaccine_strain, vaccine_tip in vaccine_tips.items()}
            tip_data[a] = summed_attribute
        distance_to_vaccine.append(tip_data)

    return distance_to_vaccine


def HI_export(process):
    prefix = os.path.join(process.config["output"]["auspice"], process.info["prefix"])

    if hasattr(process, 'HI_tree'):
        # export the raw titers
        hi_data = process.HI_tree.compile_titers()
        write_json(hi_data, prefix+'_titers.json')
        # export the tree model (avidities and potencies only)
        tree_model = {'potency':process.HI_tree.compile_potencies(),
                      'avidity':process.HI_tree.compile_virus_effects(),
                      'dTiter':{n.clade:n.dTiter for n in process.tree.tree.find_clades() if n.dTiter>1e-6}}
        write_json(tree_model, prefix+'_titer_tree_model.json')
    else:
        print('Tree model not yet trained')

    if hasattr(process, 'HI_subs'):
        # export the substitution model
        subs_model = {'potency':process.HI_subs.compile_potencies(),
                      'avidity':process.HI_subs.compile_virus_effects(),
                      'substitution':process.HI_subs.compile_substitution_effects()}
        write_json(subs_model, prefix+'_titer_subs_model.json')
    else:
        print('Substitution model not yet trained')
