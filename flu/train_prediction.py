from __future__ import division, print_function
from flu_prediction import flu_predictor
import glob
import numpy as np
from H3N2 import flu_process, H3N2_scores
from collections import defaultdict


class fitness_model_train(object):
    """docstring for fitness_model_train"""
    def __init__(self):
        super(fitness_model_train, self).__init__()

    def setup(self, prefix, years, seq_fname):
        self.pred_models = {}
        self.flu_trees = {}
        for year in years:
            time_interval = (str(year-3)+'-01-01', str(year+3)+'-01-01')
            out_specs = {'data_dir':'data/', 'prefix':'H3N2_',
                         'qualifier': time_interval[0]+'_'+time_interval[1]+'_'}
            flu = flu_process(method='SLSQP', dtps=2.0, stiffness=20, inertia=0.9,
                              time_interval=time_interval,
                              fname = seq_fname, out_specs=out_specs)
            flu.load()
            H3N2_scores(flu.tree.tree)
            flu.estimate_tree_frequencies()
            flu_pred = flu_predictor(tree=flu.tree.tree, seqs=flu.seqs,
                                     clade_frequencies = flu.tree_frequencies, years=[year],
                                     pivots=flu.tree_pivots,inertia=1.0, stiffness=20)

            flu_pred.calculate_LBI(dt=1) #dt is the max age of sequences used
            flu.raw_seqs = {}
            self.flu_trees[year] = flu
            self.pred_models[year] = flu_pred

    def add_training_frequencies(self):
        for flu_pred in self.pred_models.values():
            flu_pred.calculate_training_frequencies()

    def add_HI(self, titer_fname):
        for flu in self.flu_trees.values():
            flu.HI_model(titer_fname)

    def fit_quality(self, model, **kwarks):
        current_cost=0
        for year in self.pred_models:
            self.pred_models[year].general_fitness(model)
            current_cost += self.pred_models[year].score_model(**kwarks)
        return current_cost

    def train_model(self,model=None, method='SLSQP',**kwarks):
        '''
        train the model by optimizing the relative contributions of predictors to fitness
        arguments:
            clade_dt:   the time prior to the end of the training interval used to define clades
            model:      dictionary with predictors and intitial coefficients
            method:     method used to optimize the parameters
            metric:     metric used to calculate deviations of predictions and actual
                        frequencies. Choices are KL, sq, abs
        '''
        if model is None:
            model = {'slope':1.0}

        predictors = model.keys()
        verbose=1
        coeff = np.array([model[pred] for pred in predictors])
        from scipy.optimize import minimize

        # functon to minimize
        def cost(x):
            model = {pred:c for pred,c in zip(predictors, x)}
            current_cost = self.fit_quality(model,**kwarks)
            if verbose: print([pred+': '+str(c) for pred,c in model.iteritems()]+['dev:',current_cost])
            return current_cost

        self.sol = minimize(cost, coeff, method=method)
        self.fitness_params = {pred:c for pred,c in zip(predictors, self.sol['x'])}


    def plot_model_predictions(self, clade_dt=2, model=None):
        if model is None:
            model=self.fitness_params
        # validate by plotting trajectories
        for year in self.pred_models:
            self.pred_models[year].general_fitness(model)
            for tint in self.pred_models[year].train_intervals:
                self.pred_models[year].frequency_prediction(tint[1]-clade_dt, tint)
                self.pred_models[year].plot_prediction()
                plt.title('year: '+str(year))



if __name__=="__main__":
    import matplotlib.pyplot as plt
    plt.ion()
    #seq_fname = sorted(glob.glob('../nextstrain-db/data/flu_h3n2*fasta'))[-1]
    seq_fname = '../../nextflu2/data/flu_h3n2_gisaid.fasta'
    titer_fname = sorted(glob.glob('../nextstrain-db/data/h3n2*text'))[-1]

    fitness_trainer = fitness_model_train()
    fitness_trainer.setup('H3N2', range(1996,2012), seq_fname)
    fitness_trainer.add_training_frequencies()
    #fitness_trainer.add_HI(titer_fname)
    dt = 1.5
    fitness_trainer.train_model(model={'slope':0.4, 'LBI':0.01}, metric='abs', horizon=2, clade_dt=dt)

    fitness_trainer.plot_model_predictions(model={'LBI':0.0}, clade_dt=dt)
    fitness_trainer.plot_model_predictions(clade_dt=dt)
    for year,pred in fitness_trainer.pred_models.iteritems():
        print(year,pred.score_model(metric='abs', horizon=2, clade_dt=dt))
