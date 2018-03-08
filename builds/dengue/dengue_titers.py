def tree_additivity_symmetry(titer_model):
    '''
    The titer model makes two major assumptions:
    1 - Once we correct for virus avidity and serum potency, titers are roughly symmetric
    2 - Antigenic relationships evolve in a tree-like fashion

    We can check the validity of these assumptions by plotting:
    1 - titer symmetry (T(serum i,virus j) - vice versa)
    2 - For random quartets of antigenic distances, how do the relative distances between the tips relate to one another?
        If they're tree-like, we expect the largest - second largest distance to be roughly exponentially distributed.
        How does this compare to quartets of values pulled randomly from a normal distribution?

    Code adapted from https://github.com/blab/nextflu/blob/pnas-hi-titers/augur/src/diagnostic_figures.py#L314
    '''
    import numpy as np
    from matplotlib import pyplot as plt

    reciprocal_measurements = []
    reciprocal_measurements_titers = []
    for (testvir, serum) in titer_model.titers.titers_normalized:
        tmp_recip = [v for v in titer_model.titers.titers_normalized if serum[0]==v[0] and testvir==v[1][0]]
        for v in tmp_recip:
            val_fwd = titer_model.titers.titers_normalized[(testvir,serum)]
            val_bwd = titer_model.titers.titers_normalized[v]
            date_fwd = titer_model.node_lookup[testvir].attr['num_date']
            date_bwd = titer_model.node_lookup[serum[0]].attr['num_date']
            diff_uncorrected = val_fwd - val_bwd
            diff_corrected = (val_fwd - titer_model.serum_potency[serum] - titer_model.virus_effect[testvir])\
                            -(val_bwd - titer_model.serum_potency[v[1]] - titer_model.virus_effect[serum[0]])
            val_bwd = titer_model.titers.titers_normalized[v]
            reciprocal_measurements.append([testvir, serum, diff_uncorrected, diff_corrected, np.sign(date_fwd-date_bwd)])
            reciprocal_measurements_titers.append([testvir, serum, val_fwd, val_bwd,
                                                  (val_fwd - titer_model.serum_potency[serum] - titer_model.virus_effect[testvir]),
                                                  (val_bwd - titer_model.serum_potency[v[1]] - titer_model.virus_effect[serum[0]]),
                                                  ])
    plt.figure(figsize=(9,6))
    ax = plt.subplot(121)
    # multiple the difference by the +/- one to polarize all comparisons by date
    vals= [x[2]*x[-1] for x in reciprocal_measurements]
    plt.hist(vals, alpha=0.7, label=r"Raw $T_{ij}-T_{ji}$", normed=True)
    print("raw reciprocal titers: ", str(np.round(np.mean(vals),3))+'+/-'+str(np.round(np.std(vals),3)))
    vals= [x[3]*x[-1] for x in reciprocal_measurements]
    plt.hist(vals, alpha=0.7, label=r"Normalized $T_{ij}-T_{ji}$", normed=True)
    print("normalized reciprocal titers: ", str(np.round(np.mean(vals),3))+'+/-'+str(np.round(np.std(vals),3)))
    plt.xlabel('Titer asymmetry', fontsize=12)
    ax.tick_params(axis='both', labelsize=12)
    plt.legend(fontsize=12, handlelength=0.8)
    plt.tight_layout()

    ####  Analyze all cliques #######################################################
    all_reciprocal = list(set([v[1] for v in reciprocal_measurements_titers]))

    import networkx as nx
    from random import sample
    G = nx.Graph()
    G.add_nodes_from(all_reciprocal)
    for vi,v in enumerate(all_reciprocal):
        for w in all_reciprocal[:vi]:
            if ((v[0], w) in titer_model.titers.titers_normalized) and ((w[0], v) in titer_model.titers.titers_normalized):
                G.add_edge(v,w)
    C = nx.find_cliques(G)
    def symm_distance(v,w):
        res =  titer_model.titers.titers_normalized[(v[0], w)] - titer_model.virus_effect[v[0]] - titer_model.serum_potency[w]
        res += titer_model.titers.titers_normalized[(w[0], v)] - titer_model.virus_effect[w[0]] - titer_model.serum_potency[v]
        return res*0.5

    additivity_test = {'test':[], 'control':[]}
    n_quartets = 1000
    for clique in C:
        if len(clique)>8:
            for i in xrange(n_quartets):
                Q = sample(clique, 4)
                dists = []
                for (a,b) in [((0,1), (2,3)),((0,2), (1,3)), ((0,3), (1,2))]:
                    dists.append(symm_distance(Q[a[0]], Q[a[1]])+symm_distance(Q[b[0]], Q[b[1]]))
                dists.sort(reverse=True)
                additivity_test['test'].append(dists[0]-dists[1])

                dists = []
                for di in range(3):
                    a,b,c,d = sample(clique,4)
                    dists.append(symm_distance(a, b)+symm_distance(c,d))
                dists.sort(reverse=True)
                additivity_test['control'].append(dists[0]-dists[1])

    ax=plt.subplot(122)
    plt.hist(additivity_test['control'], alpha=0.7,normed=True, bins = np.linspace(0,3,18),
             label = 'Control, mean='+str(np.round(np.mean(additivity_test['control']),2)))
    plt.hist(additivity_test['test'], alpha=0.7,normed=True, bins = np.linspace(0,3,18),
             label = 'Quartet, mean='+str(np.round(np.mean(additivity_test['test']),2)))
    ax.tick_params(axis='both', labelsize=12)
    plt.xlabel(r'$\Delta$ top two distance sums', fontsize = 12)
    plt.legend(fontsize=12, handlelength=0.8)
    plt.tight_layout()
    plt.savefig('./processed/titer_asymmetry.png')

def titer_model(process, sanofi_strain = None, **kwargs):
    '''
    estimate a titer tree model using titers in titer_fname.
    '''
    from base.titer_model import TreeModel, SubstitutionModel
    ## TREE MODEL
    process.titer_tree = TreeModel(process.tree.tree, process.titers, **kwargs)
    process.titer_tree.prepare(**kwargs) # make training set, find subtree with titer measurements, and make_treegraph
    process.titer_tree.train(**kwargs)   # pick longest branch on path between each (test, ref) pair, assign titer drops to this branch
                                         # then calculate a cumulative antigenic evolution score for each node

    # add attributes for the estimated branch-specific titer drop values (dTiter)
    # and cumulative (from root) titer drop values (cTiter) to each branch
    for n in process.tree.tree.find_clades():
        n.attr['cTiter'] = n.cTiter
        n.attr['dTiter'] = n.dTiter

    if kwargs['training_fraction'] != 1.0:
        process.titer_tree.validate(kwargs)

    if sanofi_strain: # calculate antigenic distance from vaccine strain for each serotype-specific build
        # find the vaccine strain in the tree
        sanofi_tip = [i for i in process.tree.tree.find_clades() if i.name==sanofi_strain][0]
        # sum up dTiter on all the branches on the path between the vaccine strain and each other strain
        for tip in process.tree.tree.find_clades():
            if tip == sanofi_tip:
                tip.attr['dTiter_sanofi']=0.00
            else:
                trace = process.tree.tree.trace(tip, sanofi_tip)
                trace_dTiter = sum([i.dTiter for i in trace])
                tip.attr['dTiter_sanofi']= round(trace_dTiter, 2)

        # export for auspice visualization
        process.config["auspice"]["color_options"]["dTiter_sanofi"] = {
        "menuItem": "antigenic dist. from vaccine", "type": "continuous", "legendTitle": "log2 titer distance from sanofi vaccine strain",
        "key": "vaccine_dTiter", "vmin": "0.0", "vmax": "2.0"}

    else: # export cumulative distance for auspice vis. of the all-serotype build
        process.config["auspice"]["color_options"]["cTiter"] = {
        "menuItem": "antigenic dist. from root", "type": "continuous", "legendTitle": "log2 titer distance from root",
        "key": "cTiter", "vmin": "0.0", "vmax": "2.0"}

    if kwargs['plot_symmetry'] == True:
        tree_additivity_symmetry(process.titer_tree)

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

        validation_values = {'actual': [v[0] for v in process.titer_tree.validation.values()],
                             'predicted': [v[1] for v in process.titer_tree.validation.values()]}
        import pandas as pd
        validation_values = pd.DataFrame(validation_values)
        validation_values.to_csv(prefix+'tree_model_validation.csv')

    else:
        print('Tree model not yet trained')
