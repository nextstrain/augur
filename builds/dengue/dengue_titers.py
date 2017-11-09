def tree_additivity_symmetry(mytitermodel, mtype='tree'):
    ''' Adapted from https://github.com/blab/nextflu/blob/pnas-hi-titers/augur/src/diagnostic_figures.py#L314 '''
    import numpy as np
    from matplotlib import pyplot as plt

    figheight=8
    fs = 12
    reciprocal_measurements = []
    reciprocal_measurements_titers = []
    print(mytitermodel.__dict__.keys())
    for (testvir, serum) in mytitermodel.titers_normalized:
        tmp_recip = [v for v in mytitermodel.titers_normalized if serum[0]==v[0] and testvir==v[1][0]]
        for v in tmp_recip:
            val_fwd = mytitermodel.titers_normalized[(testvir,serum)]
            val_bwd = mytitermodel.titers_normalized[v]
            date_fwd = mytitermodel.node_lookup[testvir].attr['num_date']
            date_bwd = mytitermodel.node_lookup[serum[0]].attr['num_date']
            diff_uncorrected = val_fwd - val_bwd
            diff_corrected = (val_fwd - mytitermodel.serum_potency[serum] - mytitermodel.virus_effect[testvir])\
                            -(val_bwd - mytitermodel.serum_potency[v[1]] - mytitermodel.virus_effect[serum[0]])
            val_bwd = mytitermodel.titers_normalized[v]
            reciprocal_measurements.append([testvir, serum, diff_uncorrected, diff_corrected, np.sign(date_fwd-date_bwd)])
            reciprocal_measurements_titers.append([testvir, serum, val_fwd, val_bwd,
                                                  (val_fwd - mytitermodel.serum_potency[serum] - mytitermodel.virus_effect[testvir]),
                                                  (val_bwd - mytitermodel.serum_potency[v[1]] - mytitermodel.virus_effect[serum[0]]),
                                                  ])
    plt.figure(figsize=(1.6*figheight, figheight))
    ax = plt.subplot(121)
    plt.text(0.05, 0.93,  ('tree model' if mtype=='tree' else 'mutation model'),
             weight='bold', fontsize=fs, transform=plt.gca().transAxes)
    # multiple the difference by the +/- one to polarize all comparisons by date
    vals= [x[2]*x[-1] for x in reciprocal_measurements]
    plt.hist(vals, alpha=0.7, label=r"$H_{a\beta}-H_{b\alpha}$", normed=True)
    print("raw reciprocal titers: ", str(np.round(np.mean(vals),3))+'+/-'+str(np.round(np.std(vals),3)))
    vals= [x[3]*x[-1] for x in reciprocal_measurements]
    plt.hist(vals, alpha=0.7, label=r"sym. part", normed=True)
    print("symmetric component: ", str(np.round(np.mean(vals),3))+'+/-'+str(np.round(np.std(vals),3)))
    plt.xlabel('titer asymmetry', fontsize=fs)
    ax.tick_params(axis='both', labelsize=fs)
    # add_panel_label(ax, 'A', x_offset=-0.2)
    plt.legend(fontsize=fs, handlelength=0.8)
    plt.tight_layout()

    ####  Analyze all cliques #######################################################
    all_reciprocal = list(set([v[1] for v in reciprocal_measurements_titers]))

    import networkx as nx
    from random import sample
    G = nx.Graph()
    G.add_nodes_from(all_reciprocal)
    for vi,v in enumerate(all_reciprocal):
        for w in all_reciprocal[:vi]:
            if ((v[0], w) in mytitermodel.titers_normalized) and ((w[0], v) in mytitermodel.titers_normalized):
                G.add_edge(v,w)
    # print "generated graph of all cliques"
    C = nx.find_cliques(G)
    # print "found cliques"
    def symm_distance(v,w):
        res =  mytitermodel.titers_normalized[(v[0], w)] - mytitermodel.virus_effect[v[0]] - mytitermodel.serum_potency[w]
        res += mytitermodel.titers_normalized[(w[0], v)] - mytitermodel.virus_effect[w[0]] - mytitermodel.serum_potency[v]
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
             label = 'control, mean='+str(np.round(np.mean(additivity_test['control']),2)))
    plt.hist(additivity_test['test'], alpha=0.7,normed=True, bins = np.linspace(0,3,18),
             label = 'quartet, mean='+str(np.round(np.mean(additivity_test['test']),2)))
    ax.tick_params(axis='both', labelsize=fs)
    plt.xlabel(r'$\Delta$ top two distance sums', fontsize = fs)
    plt.legend(fontsize=fs, handlelength=0.8)
    # add_panel_label(ax, 'B', x_offset=-0.2)
    plt.tight_layout()
    plt.savefig('./processed/titer_asymmetry.png')

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

    if 'plot_symmetry' in kwargs:
        tree_additivity_symmetry(process.titer_tree)

    # SUBSTITUTION MODEL
    # process.titer_subs = SubstitutionModel(process.tree.tree, process.titers, **kwargs)
    # process.titer_subs.prepare(**kwargs)
    # process.titer_subs.train(**kwargs)

    if kwargs['training_fraction'] != 1.0:
        process.titer_tree.validate(kwargs) #(plot=True, fname='treeModel_%s.png'%lineage)
        # process.titer_subs.validate() #(plot=True, fname='subModel_%s.png'%lineage)

    process.config["auspice"]["color_options"]["cTiter"] = {
    "menuItem": "antigenic advance", "type": "continuous", "legendTitle": "log2 titer distance from root", "key": "cTiter", "vmin": "0.0", "vmax": "3.5"
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

    # if hasattr(process, 'titer_tree'):
    #     # export the substitution model
    #     titer_subs_model = {'potency':process.titer_subs.compile_potencies(),
    #                   'avidity':process.titer_subs.compile_virus_effects(),
    #                   'substitution':process.titer_subs.compile_substitution_effects()}
    #     write_json(titer_subs_model, prefix+'titer_subs_model.json')
    # else:
    #     print('Substitution model not yet trained')
