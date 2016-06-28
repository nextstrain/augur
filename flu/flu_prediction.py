from __future__ import division, print_function
import numpy as np
#from flu import H3N2
from base.prediction import tree_predictor, LBI
from H3N2 import flu_process, H3N2_scores
from collections import defaultdict

class flu_predictor(tree_predictor):
    """class to train and evaluate influenza prediction models"""
    def __init__(self, clade_frequencies=None, pivots=None, years=None, **kwarks):
        super(flu_predictor, self).__init__(**kwarks)
        self.global_freqs = clade_frequencies
        self.global_pivots = pivots
        self.eps=0.01

        # memorize all times of all nodes: necessary to determine which branches are
        # competing at a certain point in time.
        self.all_times = sorted([node.numdate for node in self.tree.find_clades()])

        # training intervals. currently use 2+2/12y worth of data up until end of Feb
        # the end of each interval will be considered the point from which on we want
        # to predict into the future
        if years is None:
            years= range(2009, 2016)
        self.train_intervals = [(x-2, x+2.0/12) for x in years]

    def calculate_training_frequencies(self, **kwarks):
        '''
        loops over all training intervals and calculates the frequency trajectories of
        clades using only time points within these intervals. the kwarks are passed
        on to the frequency estimation and could specify things like stiffness etc
        '''
        self.train_frequencies={}
        for tint in self.train_intervals:
            self.set_train(tint)
            self.train_frequencies[tint] = self.estimate_training_frequencies(**kwarks)


    def frequency_prediction(self, t_cut, train_interval):
        '''
        extrapolate clade frequencies into the future using the current fitness model
        arguments:
            t_cut:  time point used to define competing clades
            train_interval: the interval from which the future is predicted
        '''
        self.predicted_frequencies = {}
        self.t_cut=t_cut
        self.current_prediction_interval = train_interval

        t_lower = self.all_times[0]
        self.predictions = defaultdict(dict)
        if t_cut>=train_interval[1]:
            print("t_cut is past the training interval")
            return
        else:
            cut = []
            # extrapolation time point -- frequencies match exactly here.
            t0 = self.train_frequencies[train_interval][0][-2]
            # gather nodes whose branch is cut when cutting at t_cut
            tmp_freqs = {}
            for node in self.tree.find_clades():
                f0 = self.train_frequencies[train_interval][1][node.clade][-2]
                if (node.numdate>t_cut) and (node.up.numdate<t_cut):
                    cut.append(node)
                    tmp_freqs[node] = np.exp((self.global_pivots-t0)*node.fitness[train_interval])*f0

            # normalize frequencies of competing clades
            total_freq = np.sum(tmp_freqs.values(), axis=0)
            for node in cut:
                self.predictions[node] = tmp_freqs[node]/total_freq

    def general_fitness(self, model):
        '''
        calculate the fitness score of each node for each training interval
        based on the fitness model and precomputed LBI and slopes
        arguments:
            model: dictionary linking predictors to their coefficients
        '''
        for node in self.tree.find_clades():
            node.fitness={}

        for train_interval,(pivots, freq) in self.train_frequencies.iteritems():
            dt = pivots[-1]-pivots[-2]
            for node in self.tree.find_clades():
                tmp_fit = 0
                for pred, c in model.iteritems():
                    if pred == 'LBI':
                        tmp_fit += c * node.LBI[train_interval]*200
                    elif pred =='slope':
                        slope = 2.0*((freq[node.clade][-1] - freq[node.clade][-2])/dt)/\
                                    (self.eps+freq[node.clade][-1] + freq[node.clade][-2])
                        tmp_fit += c*slope
                    else:
                        tmp_fit += c*node.__getattribute__(pred)
                node.fitness[train_interval] = tmp_fit


    def train_model(self,clade_dt = 2, model=None, method='SLSQP', metric='sq', plot=False):
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
        self.calculate_LBI(dt=1) #dt is the max age of sequences used

        # functon to minimize
        def cost(x):
            self.general_fitness({pred:c for pred,c in zip(predictors, x)})
            current_cost = self.score_model(metric=metric)
            if verbose: print([pred+': '+str(c) for pred,c in zip(predictors, x)]+['dev:',current_cost])
            return current_cost

        self.sol = minimize(cost, coeff, method=method)
        self.fitness_params = {pred:c for pred,c in zip(predictors, self.sol['x'])}

        # validate by plotting trajectories
        if plot:
            self.general_fitness(self.fitness_params)
            for tint in self.train_intervals[:-1]:
                self.frequency_prediction(tint[1]-clade_dt, tint)
                self.plot_prediction()


    def prediction_error(self, metric='sq'):
        '''
        calculate the prediction error for all caldes for which predictions have been made.
        returns a vector with deviations at every time point
        '''
        train_pivots = self.train_frequencies[self.current_prediction_interval][0]
        real_freqs, pred_freqs = [], []
        for node in self.predictions:
            real_freqs.append(self.global_freqs[node.clade])
            pred_freqs.append(self.predictions[node])
        real_freqs, pred_freqs = np.array(real_freqs), np.array(pred_freqs)

        if metric=='KL': # Kullback leibler divergence
            dev = np.sum(real_freqs * np.log((self.eps+real_freqs)/(self.eps+pred_freqs)), axis=0)
        elif metric=='sq': # squared deviation
            dev = np.sum((real_freqs-pred_freqs)**2, axis=0)
        else:   #absolute deviation
            dev = np.sum(np.abs(real_freqs-pred_freqs), axis=0)

        return dev


    def score_model(self, clade_dt = 2, metric='sq', horizon=2):
        '''
        calculate the quality of the model by summing the relevant error over all
        training intervals.
        '''
        prediction_error = 0
        for tint in self.train_intervals:
            self.frequency_prediction(tint[1]-clade_dt, tint)
            # calculate deviations, restrict to future time points ((tint, tint+horizon))
            dev = self.prediction_error(metric=metric)[(self.global_pivots>tint[1])&
                                                       (self.global_pivots<tint[1]+horizon)]
            prediction_error+= np.mean(dev)
        return prediction_error


    def plot_prediction(self):
        '''
        plots the global frequencies, the predicted frequencies, and the frequencies
        in the short interval used for learning.
        '''
        from matplotlib import pyplot as plt
        import seaborn as sns
        fig, axs = plt.subplots(1,2, figsize=(12,6))

        axs[0].plot(self.t_cut*np.ones(2), [0,1], lw=3, alpha=0.3, c='k', ls='--')
        axs[0].plot(self.current_prediction_interval[1]*np.ones(2), [0,1], lw=3, alpha=0.3, c='k')

        train_pivots = self.train_frequencies[self.current_prediction_interval][0]
        train_freqs = self.train_frequencies[self.current_prediction_interval][1]
        cols = sns.color_palette()
        for node in self.predictions:
            if np.max(self.predictions[node][self.global_pivots>train_pivots[0]])>0.02:
                #print(self.predictions[t_cut_val][node])
                future_pivots = self.global_pivots>train_pivots[0]
                axs[0].plot(self.global_pivots[future_pivots],
                            self.predictions[node][future_pivots], ls='--', c=cols[node.clade%6])
                axs[0].plot(self.global_pivots, self.global_freqs[node.clade], ls='-', c=cols[node.clade%6])
                axs[0].plot(train_pivots, train_freqs[node.clade], ls='-.', c=cols[node.clade%6])

        axs[0].set_xlim(train_pivots[0]-2, train_pivots[-1]+2)
        dev = self.prediction_error()
        axs[1].plot(self.global_pivots, dev)
        axs[1].set_xlim(train_pivots[0], train_pivots[-1]+2)



    def plot_all(self):
        from matplotlib import pyplot as plt
        import seaborn as sns
        cols = sns.color_palette(n_colors=6)
        p = self.global_pivots
        fig = plt.figure(figsize = (20,7))
        ax = plt.subplot(111)
        for clade in  self.global_freqs.keys():
            f = self.global_freqs[clade]
            if np.max(f)>0.2:
                ax.plot(p, f, c=cols[clade%len(cols)], alpha=0.3, ls='--')

        for tint, (p,freq) in self.train_frequencies.iteritems():
            for clade in  freq.keys():
                if np.max(freq[clade])>0.2:
                    ax.plot(p, freq[clade], c=cols[clade%len(cols)], alpha=0.5)

    def delta_freq_vs_success(self):
        from matplotlib import pyplot as plt
        import seaborn as sns
        max_freq = {clade:(f[-1], np.max(f), np.max(np.diff(f)))
                    for clade, f in self.global_freqs.iteritems() if np.max(f)>0.03}

        dfreq_array = np.array(max_freq.values())
        plt.figure()
        plt.scatter(dfreq_array[:,0], dfreq_array[:,1])
        plt.ylabel('peak')
        plt.xlabel('final')

        plt.figure()
        plt.scatter(dfreq_array[:,0], dfreq_array[:,2])
        plt.ylabel('max jump')
        plt.xlabel('final')

        max_freq_intervals = {}

        for tint, (p,freq) in self.train_frequencies.iteritems():
            for clade in  freq.keys():
                if np.max(freq[clade])>0.1:
                    max_freq_intervals[clade] = [max_freq[clade][0], np.max(freq[clade]), np.max(np.diff(freq[clade]))]

        dfreq_array_interval = np.array(max_freq_intervals.values())
        plt.figure()
        plt.scatter(dfreq_array_interval[:,0], dfreq_array_interval[:,1])
        plt.ylabel('peak')
        plt.xlabel('final')

        plt.figure()
        plt.scatter(dfreq_array_interval[:,0], dfreq_array_interval[:,2])
        plt.ylabel('max jump')
        plt.xlabel('final')

        return dfreq_array, dfreq_array_interval




if __name__=="__main__":
    import argparse, glob
    from matplotlib import pyplot as plt
    plt.ion()

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('--load', action='store_true', help = 'recover from file')

    params = parser.parse_args()
    fname = sorted(glob.glob('../nextstrain-db/data/flu_h3n2*fasta'))[-1]
    titer_fname = sorted(glob.glob('../nextstrain-db/data/h3n2*text'))[-1]

    flu = flu_process(method='SLSQP', dtps=2.0, stiffness=20, inertia=0.9,
                      fname = fname)
    if params.load:
        flu.load()
        H3N2_scores(flu.tree.tree)
        flu.estimate_tree_frequencies()
    else:
        pass

    flu.HI_model(titer_fname)

    flu_pred = flu_predictor(tree=flu.tree.tree, seqs=flu.seqs,
                             clade_frequencies = flu.tree_frequencies,
                             pivots=flu.tree_pivots,inertia=1.0, stiffness=20)

    flu_pred.calculate_training_frequencies()


    # flu_pred.plot_all()
    # global_freqs, interval_freqs = flu_pred.delta_freq_vs_success()
#
    # failed = global_freqs[:,0]<0.1
    # fixed = global_freqs[:,0]>0.9
    # failed_peaks = global_freqs[failed,1]
    # plt.figure()
    # plt.plot(sorted(failed_peaks), np.linspace(1,0,len(failed_peaks)), label='failed')
    # plt.plot(sorted(global_freqs[:,1]), np.linspace(1,0,len(global_freqs)), label='all')
    # plt.xlabel('peak frequency')
    # plt.ylabel('fraction with higher peak')
#
    # failed_jumps = global_freqs[failed,2]
    # fixed_jumps = global_freqs[fixed,2]
    # plt.figure()
    # plt.plot(sorted(failed_jumps), np.linspace(1,0,len(failed_jumps)), label='failed')
    # plt.plot(sorted(fixed_jumps), np.linspace(1,0,len(fixed_jumps)), label='fixed')
    # plt.xlabel('max jump')
    # plt.ylabel('fraction with larger jump')
#

    flu_pred.general_fitness(flu_pred.fitness_params)
    clade_dt=2
    for tint in flu_pred.train_intervals[:-1]:
        flu_pred.frequency_prediction(tint[1]-clade_dt, tint)
        flu_pred.plot_prediction()


