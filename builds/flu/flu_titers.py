import os


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
    process.config["auspice"]["color_options"]["cTiter"] = {
        "menuItem": "antigenic advance", "type": "continuous", "legendTitle": "Antigenic Advance", "key": "cTiter"
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


def HI_export(process):
    from base.io_util import write_json
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
