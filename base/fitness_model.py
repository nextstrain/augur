import argparse
from collections import defaultdict
import numpy as np
import os
import pandas as pd
from scipy.interpolate import interp1d
from scipy.stats import linregress, pearsonr
import sys

from .io_util import write_json
from .scores import select_nodes_in_season
from .frequencies import logit_transform, KdeFrequencies
from .fitness_predictors import fitness_predictors

try:
    import itertools.izip as zip
except ImportError:
    pass

min_tips = 10
pc=1e-2
regularization = 1e-3
default_predictors = ['lb', 'ep', 'ne_star']


def process_predictor_args(predictors, params=None, sds=None):
    """Returns a predictor data structure for the given lists of predictors, params,
    and standard deviations.

    When no parameters or deviations are provided, the predictors are a simple
    list. When parameters and deviations are provided, the predictor are a
    dictionary indexed by predictor name with values corresponding to each
    predictor's param and global standard deviation.

    >>> process_predictor_args(None, None, None)
    >>> process_predictor_args(['ep'])
    ['ep']
    >>> process_predictor_args(['ep'], None, None)
    ['ep']
    >>> process_predictor_args(['ep'], [1], [5])
    {'ep': [1, 5]}
    """
    if predictors is None:
        processed_predictors = None
    elif params is None or sds is None:
        processed_predictors = predictors
    else:
        merged_params = map(list, zip(params, sds))
        processed_predictors = dict(zip(predictors, merged_params))

    return processed_predictors


def make_pivots(start, stop, pivots_per_year=12, precision=2):
    """Makes an array of pivots (i.e., timepoints) between the given start and stop
    by the given pivots per year. The generated pivots are floating point values
    that are then rounded to the requested decimal precision.

    >>> list(make_pivots(2000.0, 2001.0, 5))
    [2000.0, 2000.25, 2000.5, 2000.75, 2001.0]
    """
    # Calculate number of pivots (i.e., months) in the requested interval.
    number_of_pivots = np.ceil((stop - start) * pivots_per_year)

    # Build an evenly-spaced closed interval (including the start and stop
    # points) based on the calculated number of pivots.
    return np.around(
        np.linspace(start, stop, number_of_pivots),
        precision
    )


def matthews_correlation_coefficient(tp, tn, fp, fn):
    """Return Matthews correlation coefficient for values from a confusion matrix.
    Implementation is based on the definition from wikipedia:

    https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
    """
    numerator = (tp * tn) - (fp * fn)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denominator == 0:
            denominator = 1

    return float(numerator) / denominator


def get_matthews_correlation_coefficient_for_data_frame(freq_df):
        """Calculate Matthew's correlation coefficient from a given pandas data frame
        with columns for initial, observed, and predicted frequencies.
        """
        observed_growth = (freq_df["observed_freq"] > freq_df["initial_freq"])
        predicted_growth = (freq_df["predicted_freq"] > freq_df["initial_freq"])
        true_positives = ((observed_growth) & (predicted_growth)).sum()
        false_positives= ((~observed_growth) & (predicted_growth)).sum()

        observed_decline = (freq_df["observed_freq"] <= freq_df["initial_freq"])
        predicted_decline = (freq_df["predicted_freq"] <= freq_df["initial_freq"])
        true_negatives = ((observed_decline) & (predicted_decline)).sum()
        false_negatives = ((~observed_decline) & (predicted_decline)).sum()

        mcc = matthews_correlation_coefficient(
            true_positives,
            true_negatives,
            false_positives,
            false_negatives
        )

        return mcc


def sum_of_squared_errors(observed_freq, predicted_freq):
    """
    Calculates the sum of squared errors for observed and predicted frequencies.

    Args:
        observed_freq (numpy.ndarray): observed frequencies
        predicted_freq (numpy.ndarray): predicted frequencies

    Returns:
        float: sum of squared errors between observed and predicted frequencies
    """
    return np.sum((observed_freq - predicted_freq) ** 2)


class fitness_model(object):

    def __init__(self, tree, frequencies, predictor_input, censor_frequencies=True,
                 pivot_spacing=1.0 / 12, verbose=0, enforce_positive_predictors=True, predictor_kwargs=None,
                 cost_function=sum_of_squared_errors, delta_time=1.0, end_date=None, step_size=0.5, **kwargs):
        """

        Args:
            tree (Bio.Phylo): an annotated tree for which a fitness model is to be determined
            frequencies (KdeFrequencies): a frequency estimator and its parameters
            predictor_input: a list of predictors to fit or dict of predictors to coefficients / std deviations
            censor_frequencies (bool): whether frequencies should censor future data or not
            pivot_spacing:
            verbose:
            enforce_positive_predictors:
            predictor_kwargs:
            cost_function (callable): a function that takes observed and predicted frequencies and returns a single error value
            **kwargs:
        """
        self.tree = tree
        self.frequencies = frequencies
        self.censor_frequencies = censor_frequencies
        self.pivot_spacing = pivot_spacing
        self.verbose = verbose
        self.enforce_positive_predictors = enforce_positive_predictors
        self.estimate_coefficients = True
        self.min_freq = kwargs.get("min_freq", 0.1)
        self.max_freq = kwargs.get("max_freq", 0.99)
        self.cost_function = cost_function
        self.end_date = end_date

        if predictor_kwargs is None:
            self.predictor_kwargs = {}
        else:
            self.predictor_kwargs = predictor_kwargs

        self.time_window = kwargs.get("time_window", 6.0 / 12.0)

        if isinstance(predictor_input, dict):
            predictor_names = predictor_input.keys()
            self.estimate_coefficients = False
        else:
            predictor_names = predictor_input
        if "estimate_fitness_model" in kwargs:
            if kwargs["estimate_fitness_model"]:
                self.estimate_coefficients = True

        # Reestimate frequencies if they have not already been estimated or if internal nodes were excluded.
        if not hasattr(self.frequencies, "frequencies") or not self.frequencies.include_internal_nodes:
            sys.stderr.write("Recalculating frequencies\n")
            frequency_params = self.frequencies.get_params()
            frequency_params["include_internal_nodes"] = True
            self.frequencies = KdeFrequencies(**frequency_params)
            self.frequencies.estimate(self.tree)

        # Pivots should be defined by frequencies.
        self.pivots = self.frequencies.pivots

        # final timepoint is end of interval and is only projected forward, not tested
        self.time_interval = (self.frequencies.start_date, self.frequencies.end_date)
        self.timepoint_step_size = step_size # amount of time between timepoints chosen for fitting
        self.delta_time = delta_time         # amount of time projected forward to do fitting
        self.timepoints = np.around(
            np.append(
                make_pivots(self.time_interval[0], self.time_interval[1]-self.delta_time+0.0001, 1 / self.timepoint_step_size),
                self.time_interval[1]
            ),
            2
        )

        # Exclude the first timepoint from fitness model as censored frequencies there will be have a probability mass
        # less than 1.
        self.timepoints = self.timepoints[1:]

        # If an end date was provided, exclude all timepoints after that date.
        if self.end_date:
            original_number_of_timepoints = len(self.timepoints)
            self.timepoints = [time for time in self.timepoints
                               if time < self.end_date]
            filtered_number_of_timepoints = len(self.timepoints)
            sys.stderr.write("Fit model up to %s (filtered from %s to %s timepoints)\n" % (self.end_date, original_number_of_timepoints, filtered_number_of_timepoints))

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
                tmp_tips.append((tip_index_region_specific, node.attr["num_date"]))
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
        region = "global"
        frequencies = self.frequencies.frequencies

        # Annotate frequencies on nodes using all available data regardless of tip frequency censoring status.
        for node in self.nodes:
            node.freq = {
                region: frequencies[node.clade]
            }
            node.logit_freq = {
                region: logit_transform(frequencies[node.clade], 1e-4)
            }

        for node in self.nodes:
            interpolation = interp1d(self.rootnode.pivots, node.freq[region], kind='linear', bounds_error=True)
            node.timepoint_freqs = {}
            node.observed_final_freqs = {}
            for time in self.timepoints:
                node.timepoint_freqs[time] = np.asscalar(interpolation(time))
            for time in self.timepoints[:-1]:
                node.observed_final_freqs[time] = np.asscalar(interpolation(time + self.delta_time))

        # Estimate frequencies for tips at specific timepoints.
        # Censor future tips from estimations unless these data are explicitly allowed.
        # freq_arrays list *all* tips for each initial timepoint
        self.freq_arrays={}
        frequency_parameters = self.frequencies.get_params()
        for time in self.timepoints:
            tmp_freqs = []

            if self.censor_frequencies:
                # Censor frequencies by omitting all tips sampled after the current timepoint.
                sys.stderr.write("Calculating censored frequencies for %s\n" % time)
                frequency_parameters["max_date"] = time
                frequency_estimator = KdeFrequencies(**frequency_parameters)
                frequencies = frequency_estimator.estimate(self.tree)

                # Determine the frequency of each tip at the given timepoint.
                for tip in self.tips:
                    interpolation = interp1d(
                        self.pivots,
                        frequencies[tip.clade],
                        kind="linear",
                        bounds_error=True
                    )
                    tmp_freqs.append(np.asscalar(interpolation(time)))
            else:
                for tip in self.tips:
                    tmp_freqs.append(tip.timepoint_freqs[time])

            self.freq_arrays[time] = np.array(tmp_freqs)

    def calc_predictors(self, timepoint):
        for pred in self.predictors:
            # calculate the predictors for all nodes of the tree and save as node.attr
            if pred != 'dfreq':
                self.fp.setup_predictor(self.tree, pred, timepoint, **self.predictor_kwargs)

    def calc_time_censored_tree_frequencies(self):
        print("fitting time censored tree frequencies")
        # this doesn't interfere with the previous freq estimates via difference in region: global_censored vs global
        region = "global_censored"
        frequency_parameters = self.frequencies.get_params()
        freq_cutoff = 25.0
        pivots_fit = 6
        freq_window = 1.0
        for node in self.nodes:
            node.fit_frequencies = {}
            node.freq_slope = {}
        for time in self.timepoints:
            time_interval = [time - freq_window, time]
            node_filter_func = lambda node: node.attr['num_date'] >= time_interval[0] and node.attr['num_date'] < time_interval[1]

            # Recalculate  frequencies for the given time interval and its corresponding pivots.
            frequency_parameters["start_date"] = time_interval[0]
            frequency_parameters["end_date"] = time_interval[1]
            frequency_parameters["max_date"] = time_interval[1]
            frequency_estimator = KdeFrequencies(**frequency_parameters)
            frequencies = frequency_estimator.estimate(self.tree)

            # Annotate censored frequencies on nodes.
            # TODO: replace node-based annotation with dicts indexed by node name.
            for node in self.nodes:
                node.freq = {
                    region: frequencies[node.clade]
                }
                node.logit_freq = {
                    region: logit_transform(frequencies[node.clade], 1e-4)
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
            if self.verbose: print("calculating predictors for time", time)
            select_nodes_in_season(self.tree, time, self.time_window)
            self.calc_predictors(time)

            for node in self.nodes:
                if 'dfreq' in [x for x in self.predictors]: node.dfreq = node.freq_slope[time]

                predictors_at_time = []
                for pred in self.predictors:
                    if hasattr(node, pred):
                        predictors_at_time.append(getattr(node, pred))
                    else:
                        predictors_at_time.append(node.attr[pred])

                node.predictors[time] = np.array(predictors_at_time)

            tmp_preds = []
            for tip in self.tips:
                tmp_preds.append(tip.predictors[time])
            self.predictor_arrays[time]=np.array(tmp_preds, dtype=float)

    def standardize_predictors(self):
        self.predictor_means = {}
        self.predictor_sds = {}
        if self.verbose: print("standardizing predictors")
        for time in self.timepoints:
            values = self.predictor_arrays[time]
            weights = self.freq_arrays[time]

            # Do not calculate weighted summary statistics when all frequencies at the current timepoint are zero.
            if weights.sum() == 0:
                weights = None

            means = np.average(values, weights=weights, axis=0)
            variances = np.average((values-means)**2, weights=weights, axis=0)
            sds = np.sqrt(variances)
            self.predictor_means[time] = means
            self.predictor_sds[time] = sds

        if self.estimate_coefficients:
            self.global_sds = np.mean(self.predictor_sds.values(), axis=0)

        for time in self.timepoints:
            for node in self.nodes:
                if node.predictors[time] is not None and (self.global_sds == 0).sum() == 0:
                    node.predictors[time] = (node.predictors[time]-self.predictor_means[time]) / self.global_sds
            self.predictor_arrays[time][:,self.to_standardize] -= self.predictor_means[time][self.to_standardize]

            if (self.global_sds == 0).sum() == 0:
                self.predictor_arrays[time][:,self.to_standardize] /= self.global_sds[self.to_standardize]


    def select_clades_for_fitting(self):
        # for each time point, select clades that are within the specified frequency window
        # keep track in the dict fit_clades that maps timepoint to clade list
        self.fit_clades = {}
        for time in self.timepoints[:-1]:
            self.fit_clades[time] = []
            for node in self.nodes:
                # Only select clades for fitting if their censored frequencies are within the specified thresholds.
                node_freq = self.freq_arrays[time][node.tips].sum(axis=0)

                if self.min_freq <= node_freq <= self.max_freq:
                    # Exclude subclades whose frequency is identical to their parent clade.
                    parent_node_freq = self.freq_arrays[time][self.node_parents[node].tips].sum(axis=0)
                    if node_freq < parent_node_freq:
                        self.fit_clades[time].append(node)


    def clade_fit(self, params):
        # walk through initial/final timepoint pairs
        # tested that the sum of frequencies of tips within a clade is equal to the direct clade frequency
        timepoint_errors = []
        self.pred_vs_true = []
        pred_vs_true_values = []
        for time in self.timepoints[:-1]:
            # Project all tip frequencies forward by the specific delta time.
            all_pred_freq = self.projection(params, self.predictor_arrays[time], self.freq_arrays[time], self.delta_time)
            assert all_pred_freq.shape == self.freq_arrays[time].shape

            # normalization factor for predicted tip frequencies
            total_pred_freq = np.sum(all_pred_freq)

            # project clades forward according to strain makeup
            clade_errors = []
            tmp_pred_vs_true = []
            for clade in self.fit_clades[time]:
                # The observed final frequency is calculated for each clade from all available data.
                obs_final_freq = clade.observed_final_freqs[time]

                # The initial frequency is calculated from the sum of each clade's censored tip frequencies.
                initial_freq = self.freq_arrays[time][clade.tips].sum(axis=0)

                # The predicted final frequency is also calculated from each clade's censored tip frequencies modified
                # by the fitness and model parameters.
                pred_final_freq = np.sum(all_pred_freq[clade.tips]) / total_pred_freq

                tmp_pred_vs_true.append((initial_freq, obs_final_freq, pred_final_freq))
                pred_vs_true_values.append((time, time + self.delta_time, clade.clade, len(clade.tips), initial_freq, obs_final_freq, pred_final_freq))

            self.pred_vs_true.append(np.array(tmp_pred_vs_true))

        # Prepare a data frame with all initial, observed, and predicted frequencies by time and clade.
        self.pred_vs_true_df = pd.DataFrame(
            pred_vs_true_values,
            columns=("timepoint", "projected_timepoint", "clade", "clade_size", "initial_freq", "observed_freq", "predicted_freq")
        )

        training_error = self.cost_function(
            self.pred_vs_true_df["observed_freq"],
            self.pred_vs_true_df["predicted_freq"]
        )
        if np.isnan(training_error) or np.isinf(training_error):
            training_error = 1e10

        self.last_fit = training_error
        if self.verbose>2: print(params, self.last_fit)
        penalty = regularization*np.sum(params**2)
        if self.enforce_positive_predictors:
            for param in params:
                if param < 0:
                    penalty += 1
        return training_error + penalty

    def weighted_af(self, seqs, weights):
        af = np.zeros((4, seqs.shape[1]))
        for ni, nuc in enumerate('ACGT'):
            af[ni] += (weights*(seqs==nuc).T).sum(axis=1)/weights.sum()
        return af

    def af_fit(self, params):
        # TODO: fix me for continuous prediction
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
        if self.verbose>2: print(params, self.last_fit)
        return mean_error + regularization*np.sum(params**2)

    def fitness(self, params, pred):
        return np.sum(params*pred, axis=-1)

    def projection(self, params, pred, freqs, delta):
        return freqs * np.exp(self.fitness(params, pred) * delta);

    def minimize_clade_error(self):
        from scipy.optimize import fmin as minimizer
        if self.verbose:
            print("initial function value:", self.clade_fit(self.model_params))
            print("initial parameters:", self.model_params)
        self.model_params = minimizer(self.clade_fit, self.model_params, disp = self.verbose>1)
        if self.verbose:
            print("final function value:", self.clade_fit(self.model_params))
            print("final parameters:", self.model_params, '\n')

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
            print("initial function value:", self.af_fit(self.model_params))
            print("initial parameters:", self.model_params)
        self.model_params = minimizer(self.af_fit, self.model_params, disp = self.verbose>1)
        if self.verbose:
            print("final function value:", self.af_fit(self.model_params))
            print("final parameters:", self.model_params, '\n')


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

        print("fitting parameters of the fitness model\n")

        params_stack = []

        if self.verbose:
            print("null parameters")
        self.model_params = 0*np.ones(len(self.predictors))  # initial values
        minimize_error()
        params_stack.append((self.last_fit, self.model_params))

        for ii in xrange(niter):
            if self.verbose:
                print("iteration:", ii+1)
            self.model_params = np.random.rand(len(self.predictors)) #0*np.ones(len(self.predictors))  # initial values
            minimize_error()
            params_stack.append((self.last_fit, self.model_params))

        self.model_params = params_stack[np.argmin([x[0] for x in params_stack])][1]
        fit_func(self.model_params)
        if self.verbose:
            print("best after",niter,"iterations\nfunction value:", self.last_fit)
            print("fit parameters:")
            for pred, val in zip(self.predictors, self.model_params):
                print(pred,':', val)


    def assign_fitness(self):
        if self.verbose: print("calculating predictors for the final timepoint")
        final_timepoint = self.timepoints[-1]

        for node in self.nodes:
            if node.predictors[final_timepoint] is not None:
                node.fitness = self.fitness(self.model_params, node.predictors[final_timepoint])
            else:
                node.fitness = 0.0

            node.attr["fitness"] = node.fitness

    def assign_predicted_frequency(self, delta=1.0):
        total_freq = 0
        timepoint = self.timepoints[-1]
        for node in self.tree.get_terminals():
            pred = self.predictor_arrays[timepoint][node.tips]
            freqs = self.freq_arrays[timepoint][node.tips]
            node.predicted_freq = self.projection(self.model_params, pred, freqs, delta)[0]
            total_freq += node.predicted_freq

        for node in self.tree.get_terminals():
            node.predicted_freq /= total_freq
            node.attr["predicted_freq"] = node.predicted_freq

    def predict(self, niter = 10, estimate_frequencies = True):
        self.prep_nodes()
        self.calc_node_frequencies()
        self.calc_all_predictors(estimate_frequencies = estimate_frequencies)
        self.standardize_predictors()
        self.select_clades_for_fitting()
        if self.estimate_coefficients:
            self.learn_parameters(niter = niter, fit_func = "clade")
        self.assign_fitness()
        self.assign_predicted_frequency()

    def get_correlation(self):
        tmp = np.vstack(self.pred_vs_true)
        rho_null = pearsonr(tmp[:, 0], tmp[:, 1])
        rho_raw = pearsonr(tmp[:, 1], tmp[:, 2])
        rho_rel = pearsonr(tmp[:, 1] / tmp[:, 0],
                           tmp[:, 2] / tmp[:, 0])

        return rho_null, rho_raw, rho_rel

    def validate_prediction(self, plot=False):
        if plot:
            import matplotlib.pyplot as plt

            fig, axs = plt.subplots(1,4, figsize=(10,5))
            for time, pred_vs_true in zip(self.timepoints[:-1], self.pred_vs_true):
                # 0: initial, 1: observed, 2: predicted
                axs[0].scatter(pred_vs_true[:,1], pred_vs_true[:,2])
                axs[1].scatter(pred_vs_true[:,1]/pred_vs_true[:,0],
                               pred_vs_true[:,2]/pred_vs_true[:,0], c=pred_vs_true[0])
                for s, o, p  in pred_vs_true:
                    axs[2].arrow(s, s, o-s, p-s)
                axs[3].scatter(pred_vs_true[:,0],
                               (pred_vs_true[:,2]+0.01)/(pred_vs_true[:,1]+0.01))

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

        abs_clade_error = self.clade_fit(self.model_params)
        print("Abs clade error:", abs_clade_error)

        rho_null, rho_raw, rho_rel = self.get_correlation()
        print("Pearson's R, null:", rho_null)
        print("Pearson's R, raw:", rho_raw)
        print("Pearson's R, rel:", rho_rel)

        # pred_vs_true is initial, observed, predicted
        tmp = np.vstack(self.pred_vs_true)

        growth_list = [pred > initial for (initial, obs, pred) in tmp if obs > initial]
        correct_growth = growth_list.count(True)
        total_growth = float(len(growth_list))
        decline_list = [pred < initial for (initial, obs, pred) in tmp if obs < initial]
        correct_decline = decline_list.count(True)
        total_decline = float(len(decline_list))

        trajectory_mcc = get_matthews_correlation_coefficient_for_data_frame(self.pred_vs_true_df)

        print("Correct at predicting growth: %s (%s / %s)" % ((correct_growth / total_growth), correct_growth, total_growth))
        print("Correct at predicting decline: %s (%s / %s)" % ((correct_decline / total_decline), correct_decline, total_decline))
        print("Correct classification:",  (correct_growth+correct_decline) / (total_growth+total_decline))
        print("Matthew's correlation coefficient: %s" % trajectory_mcc)
        print("Params:")
        print(zip(self.predictors, self.model_params))

        pred_data = []
        for time, pred_vs_true in zip(self.timepoints[:-1], self.pred_vs_true):
            for entry in pred_vs_true:
                pred_data.append(np.append(entry, time))
        pred_vs_true_df = pd.DataFrame(pred_data, columns=['initial', 'obs', 'pred', 'time'])

        output_dir = "data"
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        pred_vs_true_df.to_csv(os.path.join(output_dir, "prediction_pairs.tsv"), sep="\t", index=False)

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
                interpolation = interp1d(self.rootnode.pivots, clade.freq['global'], kind='linear', bounds_error=True)
                for delta in np.arange(-0.5, 1.1, 0.1):
                    if time + delta >= self.rootnode.pivots[0] and time + delta <= self.rootnode.pivots[-1]:
                        obs_freq = np.asscalar(interpolation(time+delta))
                        pred_freq = obs_freq
                        if delta >= 0:
                            total_pred_freq = np.sum(self.projection(self.model_params, all_pred, all_freqs, delta))
                            pred_freq = np.sum(self.projection(self.model_params, pred, freqs, delta)) / total_pred_freq
                        self.trajectory_data.append([series, str(clade), time, time+delta, obs_freq, pred_freq])
                series += 1

        self.trajectory_data_df = pd.DataFrame(self.trajectory_data, columns=['series', 'clade', 'initial_time', 'time', 'obs', 'pred'])
        self.trajectory_data_df.to_csv("data/prediction_trajectories.tsv", sep="\t", index=False)

        import seaborn as sns
        import matplotlib.pyplot as plt
        cols = sns.color_palette(n_colors=6)
        fig, axs = plt.subplots(6,4, sharey=True)
        for tp, ax in zip(self.timepoints[:-1], axs.flatten()):
            traj = self.trajectory_data_df[self.trajectory_data_df.initial_time == tp]
            clades = np.unique(traj['series'])
            for ci in clades:
                tmp = traj[traj['series']==ci]
                ax.plot(tmp['time'], tmp['obs'], ls='-', c=cols[ci%6])
                ax.plot(tmp['time'], tmp['pred'], ls='--', c=cols[ci%6])

    def to_data_frame(self):
        """Return a data frame representing the fitness model's inputs for each tip.

        Inputs include tip metadata, frequencies, and standardized predictor values.
        """
        records = []

        # Include only timepoints used to fit the model itself (excluding the last timepoint).
        for timepoint in self.timepoints[:-1]:
            # Create a record for each clade fit by the model despite nesting of clades.
            for clade in self.fit_clades[timepoint]:
                # Store information for each tip in the current clade despite
                # redundancy of information. This enables refitting the model
                # with the same data later.
                for tip_index in clade.tips:
                    tip = self.tips[tip_index]
                    record = {
                        "timepoint": timepoint,
                        "clade_name": clade.name,
                        "name": tip.name,
                        "frequency": self.freq_arrays[timepoint][tip_index]
                    }

                    # Store standardized predictor values using a column name
                    # prefix that enables downstream analyses to easily identify
                    # predictor columns.
                    for predictor_index, predictor in enumerate(self.predictors):
                        record["predictor:%s" % predictor] = self.predictor_arrays[timepoint][tip_index][predictor_index]

                    records.append(record)

        return pd.DataFrame(records)

    def to_json(self, filename):
        """Export fitness model parameters, data, and accuracy statistics to JSON.
        """
        # Convert predictor parameters to a data frame to easily export as
        # records.
        params_df = pd.DataFrame({
            "predictor": self.predictors,
            "param": self.model_params.tolist(),
            "global_sd": self.global_sds.tolist()
        })

        correlation_null, correlation_raw, correlation_rel = self.get_correlation()
        mcc = get_matthews_correlation_coefficient_for_data_frame(self.pred_vs_true_df)

        # Do not try to export titer data if it was provided to the model.
        predictor_kwargs = self.predictor_kwargs.copy()
        if "titers" in predictor_kwargs:
            del predictor_kwargs["titers"]

        data = {
            "params": params_df.to_dict(orient="records"),
            "predictor_kwargs": predictor_kwargs,
            "data": self.pred_vs_true_df.to_dict(orient="records"),
            "accuracy": {
                "clade_error": self.clade_fit(self.model_params),
                "correlation_rel": correlation_rel[0],
                "mcc": mcc
            },
            "delta_time": self.delta_time,
            "step_size": self.timepoint_step_size,
            "end_date": self.end_date
        }

        predictor_arrays = {}
        for key in self.predictor_arrays:
            predictor_arrays[key] = self.predictor_arrays[key].tolist()

        data["predictor_arrays"] = predictor_arrays

        freq_arrays = {}
        for key in self.freq_arrays:
            freq_arrays[key] = self.freq_arrays[key].tolist()

        data["freq_arrays"] = freq_arrays

        write_json(data, filename)


def main(params):
    import time
    from io_util import read_json
    from io_util import write_json
    from tree_util import json_to_dendropy, dendropy_to_json

    print("--- Start fitness model optimization at " + time.strftime("%H:%M:%S") + " ---")

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
