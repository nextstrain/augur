import argparse
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import linregress
from collections import defaultdict
from itertools import izip

from frequencies import logit_transform, make_pivots, tree_frequencies
from fitness_predictors import fitness_predictors

min_freq = 0.1
max_freq = 0.99
min_tips = 10
pc=1e-2
regularization = 1e-3
default_predictors = ['lb', 'ep', 'ne_star']

class fitness_model(object):

    def __init__(self, tree, frequencies, predictor_input = ['ep', 'lb', 'dfreq'], verbose = 0, enforce_positive_predictors = True, **kwargs):
        '''
        parameters:
        tree -- tree of sequences for which a fitness model is to be determined
        frequencies -- dictionary of precalculated clade frequencies indexed by region (e.g., "global")
        predictor_input -- list of predictors to fit or dict of predictors to coefficients / std deviations
        '''
        self.tree = tree
        self.frequencies = frequencies
        self.verbose = verbose
        self.enforce_positive_predictors = enforce_positive_predictors
        self.estimate_coefficients = True

        if isinstance(predictor_input, dict):
            predictor_names = predictor_input.keys()
            self.estimate_coefficients = False
        else:
            predictor_names = predictor_input
        if "estimate_fitness_model" in kwargs:
            if kwargs["estimate_fitness_model"]:
                self.estimate_coefficients = True

        # If a time interval is provided in the kwargs, set it as an
        # attribute. Otherwise, assume it is set by a parent or mix-in
        # class.
        if not hasattr(self, "time_interval") and "time_interval" in kwargs:
            self.time_interval = kwargs["time_interval"]

        # Similarly, if pivots per year are given in kwargs, assign them
        # here.
        if not hasattr(self, "pivots_per_year") and "pivots_per_year" in kwargs:
            self.pivots_per_year = kwargs["pivots_per_year"]

        # If pivots have not been calculated yet, calculate them here.
        if (not hasattr(self, "pivots") and
            hasattr(self, "time_interval") and
            hasattr(self, "pivots_per_year")):
            self.pivots = make_pivots(self.pivots_per_year, self.time_interval)

        # final timepoint is end of interval and is only projected forward, not tested
        self.timepoint_step_size = 0.5      # amount of time between timepoints chosen for fitting
        self.delta_time = 1.0               # amount of time projected forward to do fitting
        self.timepoints = np.append(np.arange(self.time_interval[0], self.time_interval[1]-self.delta_time+0.0001, self.timepoint_step_size), self.time_interval[1])

        self.predictors = predictor_names

        self.model_params = np.zeros(len(self.predictors))
        if isinstance(predictor_input, dict):
            self.model_params = np.array([predictor_input[k][0] for k in predictor_names])

        self.to_standardize = np.array([p!='dfreq' for p in self.predictors])
        if isinstance(predictor_input, dict):
            self.global_sds = np.array([predictor_input[k][1] for k in predictor_names])
        else:
            self.global_sds = np.zeros(len(self.predictors))

        self.fp = fitness_predictors(predictor_names = predictor_names, **kwargs)

        # Map node names to parents.
        self.node_parents = {}
        for clade in self.tree.find_clades(order='level'):
            for child in clade:
                self.node_parents[child] = clade

    def prep_nodes(self):
        """Assigns data from the tree to top-level fitness model attributes.

        TODO: consider moving this code directly into the `predict`
        method since it is only ever called there.
        """
        self.nodes = [node for node in self.tree.find_clades(order="postorder")]
        self.tips = [node for node in self.nodes if node.is_terminal()]
        self.rootnode = self.tree.root
        self.rootnode.pivots = self.pivots

        # Create a list of tip indices under node.tips that map to self.tips
        # list.
        tip_index_region_specific = 0
        for node in self.nodes:
            tmp_tips = []
            if node.is_terminal():
                tmp_tips.append((tip_index_region_specific, node.numdate))
                tip_index_region_specific += 1

            for child in node.clades:
                tmp_tips.extend(child.tips)

            # Sort tips by corresponding date.
            node.tips = np.array([x for x in sorted(tmp_tips, key = lambda x: x[1])])

        # Erase the dates from the tip lists and cast to int such that they can
        # be used for indexing. These operations must happen after all nodes
        # have been processed and sorted.
        for node in self.nodes:
            if len(node.tips.shape) == 2:
                node.tips = np.array(node.tips[:, 0], dtype=int)
            else:
                node.tips = np.array([], dtype=int)

    def calc_node_frequencies(self):
        '''
        goes over all nodes and calculates frequencies at timepoints based on previously calculated frequency trajectories
        '''
        print("calculating global node frequencies")
        region = "global"

        # Calculate global tree/clade frequencies if they have not been calculated already.
        if region not in self.frequencies or self.rootnode.clade not in self.frequencies["global"]:
            tree_freqs = tree_frequencies(self.tree, self.pivots)
            tree_freqs.estimate_clade_frequencies()
            self.frequencies[region] = tree_freqs.frequencies

        # Annotate frequencies on nodes.
        # TODO: replace node-based annotation with dicts indexed by node name.
        for node in self.nodes:
            node.freq = {
                region: self.frequencies[region][node.clade]
            }
            node.logit_freq = {
                region: logit_transform(self.frequencies[region][node.clade], 1e-4)
            }

        for node in self.nodes:
            interpolation = interp1d(self.rootnode.pivots, node.freq[region], kind='linear', bounds_error=False)
            node.timepoint_freqs = defaultdict(float)
            node.delta_freqs = defaultdict(float)
            for time in self.timepoints:
                node.timepoint_freqs[time] = np.asscalar(interpolation(time))
            for time in self.timepoints[:-1]:
                node.delta_freqs[time] = np.asscalar(interpolation(time + self.delta_time))
        # freq_arrays list *all* tips for each initial timepoint
        self.freq_arrays={}
        for time in self.timepoints:
            tmp_freqs = []
            for tip in self.tips:
                tmp_freqs.append(tip.timepoint_freqs[time])
            self.freq_arrays[time] = np.array(tmp_freqs)


    def calc_predictors(self, timepoint):
        for pred in self.predictors:
            # calculate the predictors for all nodes of the tree and save as node.attr
            if pred != 'dfreq':
                self.fp.setup_predictor(self.tree, pred, timepoint)

    def select_nodes_in_season(self, timepoint):
        # used by fitness_predictors:calc_LBI and fitness_predictors:calc_epitope_cross_immunity
        # TODO: fix me for continous time model
        cutoff = 0.001
#           for node in self.nodes:
#               #if season in node.season_tips and len(node.season_tips[season])>0:
#               if node.timepoint_freqs[timepoint] > cutoff:
#                   node.alive=True
#               else:
#                   node.alive=False
        for node in self.tree.find_clades(order="postorder"):
            #if season in node.season_tips and len(node.season_tips[season])>0:
            if not node.is_terminal():
                node.alive = any(ch.alive for ch in node.clades)
            else:
                if node.numdate <= timepoint and node.numdate > timepoint - 6.0/12.0:
                    node.alive=True
                else:
                    node.alive=False

    def calc_time_censored_tree_frequencies(self):
        print("fitting time censored tree frequencies")
        # this doesn't interfere with the previous freq estimates via difference in region: global_censored vs global
        region = "global_censored"
        freq_cutoff = 25.0
        total_pivots = 12
        pivots_fit = 6
        freq_window = 1.0
        for node in self.nodes:
            node.fit_frequencies = {}
            node.freq_slope = {}
        for time in self.timepoints:
            time_interval = [time - freq_window, time]
            pivots = make_pivots(total_pivots, time_interval)
            node_filter_func = lambda node: node.numdate >= time_interval[0] and node.numdate < time_interval[1]

            # Recalculate tree frequencies for the given time interval and its
            # corresponding pivots.
            tree_freqs = tree_frequencies(self.tree, pivots, node_filter=node_filter_func)
            tree_freqs.estimate_clade_frequencies()
            self.frequencies[time] = tree_freqs.frequencies

            # Annotate censored frequencies on nodes.
            # TODO: replace node-based annotation with dicts indexed by node name.
            for node in self.nodes:
                node.freq = {
                    region: self.frequencies[time][node.clade]
                }
                node.logit_freq = {
                    region: logit_transform(self.frequencies[time][node.clade], 1e-4)
                }

            for node in self.nodes:
                if node.logit_freq[region] is not None:
                    node.fit_frequencies[time] = np.minimum(freq_cutoff, np.maximum(-freq_cutoff,node.logit_freq[region]))
                else:
                    node.fit_frequencies[time] = self.node_parents[node].fit_frequencies[time]
                try:
                    slope, intercept, rval, pval, stderr = linregress(pivots[pivots_fit:-1], node.fit_frequencies[time][pivots_fit:-1])
                    node.freq_slope[time] = slope
                except:
                    import ipdb; ipdb.set_trace()
        # reset pivots in tree to global pivots
        self.rootnode.pivots = self.pivots


    def calc_all_predictors(self, estimate_frequencies = True):
        if estimate_frequencies and 'dfreq' in [x for x in self.predictors]:
            self.calc_time_censored_tree_frequencies()
        # predictor_arrays list *all* tips for each timepoint
        self.predictor_arrays={}
        for node in self.nodes:
            node.predictors = {}
        for time in self.timepoints:
            if self.verbose: print "calculating predictors for time", time
            self.select_nodes_in_season(time)
            self.calc_predictors(time)
            for node in self.nodes:
                if 'dfreq' in [x for x in self.predictors]: node.dfreq = node.freq_slope[time]
                node.predictors[time] = np.array([node.__getattribute__(pred) for pred in self.predictors])
            tmp_preds = []
            for tip in self.tips:
                tmp_preds.append(tip.predictors[time])
            self.predictor_arrays[time]=np.array(tmp_preds, dtype=float)

    def standardize_predictors(self):
        self.predictor_means = {}
        self.predictor_sds = {}
        if self.verbose: print "standardizing predictors"
        for time in self.timepoints:
            values = self.predictor_arrays[time]
            weights = self.freq_arrays[time]
            means = np.average(values, weights=weights, axis=0)
            variances = np.average((values-means)**2, weights=weights, axis=0)
            sds = np.sqrt(variances)
            self.predictor_means[time] = means
            self.predictor_sds[time] = sds

        if self.estimate_coefficients:
            self.global_sds = np.mean(self.predictor_sds.values(), axis=0)

        for time in self.timepoints:
            for node in self.nodes:
                if node.predictors[time] is not None:
                    node.predictors[time] = (node.predictors[time]-self.predictor_means[time]) / self.global_sds
            self.predictor_arrays[time][:,self.to_standardize] -= self.predictor_means[time][self.to_standardize]
            self.predictor_arrays[time][:,self.to_standardize] /= self.global_sds[self.to_standardize]


    def select_clades_for_fitting(self):
        # for each time point, select clades that are within the specified frequency window
        # keep track in the dict fit_clades that maps timepoint to clade list
        self.fit_clades = {}
        for time in self.timepoints[:-1]:
            self.fit_clades[time] = []
            for node in self.nodes:
                if  node.timepoint_freqs[time] >= min_freq and \
                    node.timepoint_freqs[time] <= max_freq and \
                    node.timepoint_freqs[time] < self.node_parents[node].timepoint_freqs[time]:
                    self.fit_clades[time].append(node)


    def clade_fit(self, params):
        # walk through initial/final timepoint pairs
        # tested that the sum of frequencies of tips within a clade is equal to the direct clade frequency
        timepoint_errors = []
        self.pred_vs_true = []
        for time in self.timepoints[:-1]:

            # normalization factor for predicted tip frequencies
            total_pred_freq = np.sum(self.projection(params, self.predictor_arrays[time], self.freq_arrays[time], self.delta_time))

            # project clades forward according to strain makeup
            clade_errors = []
            tmp_pred_vs_true = []
            for clade in self.fit_clades[time]:
                initial_freq = clade.timepoint_freqs[time]
                obs_final_freq = clade.delta_freqs[time]
                pred = self.predictor_arrays[time][clade.tips]
                freqs = self.freq_arrays[time][clade.tips]
                pred_final_freq = np.sum(self.projection(params, pred, freqs, self.delta_time)) / total_pred_freq
                tmp_pred_vs_true.append((initial_freq, obs_final_freq, pred_final_freq))
                clade_errors.append(np.absolute(pred_final_freq - obs_final_freq))
            timepoint_errors.append(np.mean(clade_errors))
            self.pred_vs_true.append(np.array(tmp_pred_vs_true))

        mean_error = np.mean(timepoint_errors)
        if any(np.isnan(timepoint_errors)+np.isinf(timepoint_errors)):
            mean_error = 1e10
        self.last_fit = mean_error
        if self.verbose>2: print params, self.last_fit
        penalty = regularization*np.sum(params**2)
        if self.enforce_positive_predictors:
            for param in params:
                if param < 0:
                    penalty += 1
        return mean_error + penalty

    def weighted_af(self, seqs, weights):
        af = np.zeros((4, seqs.shape[1]))
        for ni, nuc in enumerate('ACGT'):
            af[ni] += (weights*(seqs==nuc).T).sum(axis=1)/weights.sum()
        return af

    def af_fit(self, params):
        # TODO: fix me for continuos prediction
        seasonal_errors = []
        self.pred_vs_true = []
        for s,t in self.fit_test_season_pairs:
            weights = np.exp(self.fitness(params, self.predictor_arrays[s][self.tree.root.season_tips[s],:]))
            pred_af = self.weighted_af(self.seqs[s],weights)
            #seasonal_errors.append(np.mean(np.sum((pred_af-self.af[t])**2, axis=0), axis=0))
            future_diameter = 0.5*np.sum(np.sum(self.af[t]*(1-self.af[t]), axis=0), axis=0)
            seasonal_errors.append(np.sum(np.sum(pred_af*(1-self.af[t]), axis=0), axis=0)-future_diameter)
            good_ind = self.af[s]*(1-self.af[s])>0.05
            self.pred_vs_true.append(np.array(zip(self.af[s][good_ind], self.af[t][good_ind], pred_af[good_ind])))

        mean_error = np.mean(seasonal_errors)
        if any(np.isnan(seasonal_errors)+np.isinf(seasonal_errors)):
            mean_error = 1e10
        self.last_fit = mean_error
        if self.verbose>2: print params, self.last_fit
        return mean_error + regularization*np.sum(params**2)

    def fitness(self, params, pred):
        return np.sum(params*pred, axis=-1)

    def projection(self, params, pred, freqs, delta):
        return freqs * np.exp(self.fitness(params, pred) * delta);

    def minimize_clade_error(self):
        from scipy.optimize import fmin as minimizer
        if self.verbose:
            print "initial function value:", self.clade_fit(self.model_params)
            print "initial parameters:", self.model_params
        self.model_params = minimizer(self.clade_fit, self.model_params, disp = self.verbose>1)
        if self.verbose:
            print "final function value:", self.clade_fit(self.model_params)
            print "final parameters:", self.model_params, '\n'

    def prep_af(self):
        if not hasattr(self,'variable_nuc'):
            self.determine_variable_positions()
        fit_aln = np.zeros((len(self.tips), len(self.variable_nuc)), dtype='S1')
        for i in range(len(self.tips)):
            tip = self.tips[i]
            fit_aln[i] = np.fromstring(tip.seq, 'S1')[self.variable_nuc]
        self.seqs = fit_aln
        self.af = {}
        for time in self.timepoints:
            self.af[time] = self.weighted_af(self.seqs, self.freq_arrays[time])

    def minimize_af_error(self):
        from scipy.optimize import fmin as minimizer
        if self.verbose:
            print "initial function value:", self.af_fit(self.model_params)
            print "initial parameters:", self.model_params
        self.model_params = minimizer(self.af_fit, self.model_params, disp = self.verbose>1)
        if self.verbose:
            print "final function value:", self.af_fit(self.model_params)
            print "final parameters:", self.model_params, '\n'


    def learn_parameters(self, niter = 10, fit_func = "clade"):
        if fit_func=='clade':
            minimize_error=self.minimize_clade_error
            fit_func=self.clade_fit
        elif fit_func=="af":
            minimize_error=self.minimize_af_error
            fit_func=self.af_fit
        else:
            print("fit function", fit_func,"does not exist")
            raise NotImplementedError

        print "fitting parameters of the fitness model\n"

        params_stack = []

        if self.verbose:
            print "null parameters"
        self.model_params = 0*np.ones(len(self.predictors))  # initial values
        minimize_error()
        params_stack.append((self.last_fit, self.model_params))

        for ii in xrange(niter):
            if self.verbose:
                print "iteration:", ii+1
            self.model_params = np.random.rand(len(self.predictors)) #0*np.ones(len(self.predictors))  # initial values
            minimize_error()
            params_stack.append((self.last_fit, self.model_params))

        self.model_params = params_stack[np.argmin([x[0] for x in params_stack])][1]
        fit_func(self.model_params)
        if self.verbose:
            print "best after",niter,"iterations\nfunction value:", self.last_fit
            print "fit parameters:"
            for pred, val in izip(self.predictors, self.model_params):
                print pred,':', val


    def assign_fitness(self):
        if self.verbose: print "calculating predictors for the final timepoint"
        final_timepoint = self.timepoints[-1]

        for node in self.nodes:
            if node.predictors[final_timepoint] is not None:
                node.fitness = self.fitness(self.model_params, node.predictors[final_timepoint])
            else:
                node.fitness = 0.0

    def predict(self, niter = 10, estimate_frequencies = True):
        self.prep_nodes()
        self.calc_node_frequencies()
        self.calc_all_predictors(estimate_frequencies = estimate_frequencies)
        self.standardize_predictors()
        self.select_clades_for_fitting()
        if self.estimate_coefficients:
            self.learn_parameters(niter = niter, fit_func = "clade")
        self.assign_fitness()

    def validate_prediction(self):
        import matplotlib.pyplot as plt
        from scipy.stats import spearmanr
        fig, axs = plt.subplots(1,4, figsize=(10,5))
        for time, pred_vs_true in izip(self.timepoints[:-1], self.pred_vs_true):
            # 0: initial, 1: observed, 2: predicted
            axs[0].scatter(pred_vs_true[:,1], pred_vs_true[:,2])
            axs[1].scatter(pred_vs_true[:,1]/pred_vs_true[:,0],
                           pred_vs_true[:,2]/pred_vs_true[:,0], c=pred_vs_true[0])
            for s, o, p  in pred_vs_true:
                axs[2].arrow(s, s, o-s, p-s)
            axs[3].scatter(pred_vs_true[:,0],
                           (pred_vs_true[:,2]+0.01)/(pred_vs_true[:,1]+0.01))

        # pred_vs_true is initial, observed, predicted
        tmp = np.vstack(self.pred_vs_true)
        print("Abs clade error:"), self.clade_fit(self.model_params)
        print("Spearman's rho, null:", spearmanr(tmp[:,0], tmp[:,1]))
        print("Spearman's rho, raw:", spearmanr(tmp[:,1], tmp[:,2]))
        print("Spearman's rho, rel:", spearmanr(tmp[:,1]/tmp[:,0],
                                              tmp[:,2]/tmp[:,0]))

        growth_list = [pred > initial for (initial, obs, pred) in tmp if obs > initial]
        correct_growth = growth_list.count(True)
        total_growth = float(len(growth_list))
        decline_list = [pred < initial for (initial, obs, pred) in tmp if obs < initial]
        correct_decline = decline_list.count(True)
        total_decline = float(len(decline_list))

        print ("Correct at predicting growth:", correct_growth / total_growth)
        print ("Correct at predicting decline:",  correct_decline / total_decline)
        print ("Correct classification:",  (correct_growth+correct_decline) / (total_growth+total_decline))

        axs[0].set_ylabel('predicted')
        axs[0].set_xlabel('observed')
        axs[1].set_ylabel('predicted/initial')
        axs[1].set_xlabel('observed/initial')
        axs[1].set_yscale('linear')
        axs[1].set_xscale('linear')
        axs[2].set_ylabel('predicted')
        axs[2].set_xlabel('observed')
        axs[2].set_ylim(-0.1, 1.1)
        axs[2].set_xlim(-0.1, 1.1)
        axs[3].set_ylabel('predicted / observed')
        axs[3].set_xlabel('initial')
        axs[3].set_yscale('log')

        import pandas as pd
        pred_data = []
        for time, pred_vs_true in izip(self.timepoints[:-1], self.pred_vs_true):
            for entry in pred_vs_true:
                pred_data.append(np.append(entry, time))
        self.pred_vs_true_df = pd.DataFrame(pred_data, columns=['initial', 'obs', 'pred', 'time'])
        self.pred_vs_true_df.to_csv("data/prediction_pairs.tsv", sep="\t", index=False)

    def validate_trajectories(self):
        '''
        Project clade trajectories based on fitted fitness model and compare to observed trajectories
        '''
        self.trajectory_data = []
        series = 0
        for time in self.timepoints[:-1]:
            all_pred = self.predictor_arrays[time]
            all_freqs = self.freq_arrays[time]
            for clade in self.fit_clades[time]:
                initial_freq = clade.timepoint_freqs[time]
                pred = all_pred[clade.tips]
                freqs = all_freqs[clade.tips]
                interpolation = interp1d(self.rootnode.pivots, clade.freq['global'], kind='linear', bounds_error=False)
                for delta in np.arange(-0.5, 1.1, 0.1):
                    if time + delta >= self.rootnode.pivots[0] and time + delta <= self.rootnode.pivots[-1]:
                        obs_freq = np.asscalar(interpolation(time+delta))
                        pred_freq = obs_freq
                        if delta >= 0:
                            total_pred_freq = np.sum(self.projection(self.model_params, all_pred, all_freqs, delta))
                            pred_freq = np.sum(self.projection(self.model_params, pred, freqs, delta)) / total_pred_freq
                        self.trajectory_data.append([series, str(clade), time, time+delta, obs_freq, pred_freq])
                series += 1

        import pandas as pd
        self.trajectory_data_df = pd.DataFrame(self.trajectory_data, columns=['series', 'clade', 'initial_time', 'time', 'obs', 'pred'])
        self.trajectory_data_df.to_csv("data/prediction_trajectories.tsv", sep="\t", index=False)
        import bokeh.charts as bk
        bk.defaults.height = 250
        bk.output_file("lines.html", title="line plot example")
        lines = []
        for time in self.timepoints[:-1]:
            line = bk.Line(self.trajectory_data_df[self.trajectory_data_df.initial_time == time],
                    x='time', y=['obs', 'pred'], dash=['obs', 'pred'], color='clade',
                    xlabel='Date', ylabel='Frequency', tools=False)
            lines.append(line)
        bk.show(bk.vplot(*lines))
#           import seaborn as sns
#           import matplotlib.pyplot as plt
#           cols = sns.color_palette(n_colors=6)
#           fig, axs = plt.subplots(6,4, sharey=True)
#           for tp, ax in zip(self.timepoints[:-1], axs.flatten()):
#               traj = self.trajectory_data_df[self.trajectory_data_df.initial_time == tp]
#               clades = np.unique(traj['series'])
#               for ci in clades:
#                   tmp = traj[traj['series']==ci]
#                   ax.plot(tmp['time'], tmp['obs'], ls='-', c=cols[ci%6])
#                   ax.plot(tmp['time'], tmp['pred'], ls='--', c=cols[ci%6])

def test(params):
    from io_util import read_json
    from tree_util import json_to_dendropy, to_Biopython, color_BioTree_by_attribute
    from Bio import Phylo
    tree_fname='data/tree_refine_10y_50v.json'
    tree =  json_to_dendropy(read_json(tree_fname))
    fm = fitness_model(tree, predictors = params['predictors'], verbose=2)
    fm.predict(niter = params['niter'])
    #btree = to_Biopython(tree)
    #color_BioTree_by_attribute(btree, 'fitness')
    #Phylo.draw(btree, label_func=lambda x:'')
    return fm

def main(params):
    import time
    from io_util import read_json
    from io_util import write_json
    from tree_util import json_to_dendropy, dendropy_to_json

    print "--- Start fitness model optimization at " + time.strftime("%H:%M:%S") + " ---"

    tree_fname='data/tree_refine.json'
    tree =  json_to_dendropy(read_json(tree_fname))
    fm = fitness_model(tree, predictors = params['predictors'], verbose=1)
    fm.predict(niter = params['niter'])
    out_fname = "data/tree_fitness.json"
    write_json(dendropy_to_json(tree.root), out_fname)
    return out_fname

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Optimize predictor coefficients')
    parser.add_argument('-n', '--niter', type = int, default=10, help='number of replicate optimizations')
    parser.add_argument("-t", "--test", help="run test", action="store_true")
    parser.add_argument('-p', '--predictors', default=default_predictors, help='predictors to optimize', nargs='+')
    params = parser.parse_args().__dict__
    if params['test']:
        fm = test(params)
    else:
        main(params)
