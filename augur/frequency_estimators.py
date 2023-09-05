# estimates clade frequencies
from __future__ import division, print_function
from collections import deque
import datetime
import isodate
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import norm
import sys

from .dates import numeric_date

debug = False
log_thres = 10.0


class TreeKdeFrequenciesError(Exception):
    """Represents an error estimating KDE frequencies for a tree.
    """
    pass


def get_pivots(observations, pivot_interval, start_date=None, end_date=None, pivot_interval_units="months"):
    """Calculate pivots for a given list of floating point observation dates and
    interval between pivots.

    Start and end pivots will be based on the range of given observed dates,
    unless a start or end date are provided to override these defaults. Pivots
    include the end date and, where possible for the given interval, the start
    date.

    Parameters
    ----------
    observations : list
        a list of observed floating point dates per sample
    pivot_interval : int
        number of months (or weeks) between pivots
    start_date : float
        optional start of the pivots interval
    end_date : float
        optional end of the pivots interval
    pivot_interval : str
        whether pivots are measured in "months" or in "weeks"

    Returns
    -------
    pivots : numpy.ndarray
        floating point pivots spanning the given the dates

    """
    # Convert months between pivots to pivot frequency.
    pivot_frequency = pivot_interval / 12.0
    if pivot_interval_units == "weeks":
        pivot_frequency = pivot_interval / 52.1429

    pivot_start = start_date if start_date else np.floor(np.min(observations) / pivot_frequency) * pivot_frequency
    pivot_end = end_date if end_date else np.ceil(np.max(observations) / pivot_frequency) * pivot_frequency

    if pivot_interval_units == "months":
        duration_str = f'P{pivot_interval}M'
    elif pivot_interval_units == "weeks":
        duration_str = f'P{pivot_interval}W'
    else:
        raise ValueError(f"The given interval unit '{pivot_interval_units}' is not supported.")

    # Construct standard Python date instances from numeric dates via ISO-8601
    # dates and the corresponding delta time for the interval between pivots.
    start = datetime.datetime.strptime(float_to_datestring(pivot_start), "%Y-%m-%d")
    end = datetime.datetime.strptime(float_to_datestring(pivot_end), "%Y-%m-%d")
    delta = isodate.parse_duration(duration_str)

    # Start calculating pivots from the end date (inclusive), working backwards
    # in time by a time delta that matches the user-requested interval. Include
    # the start date in the final pivots when the interval spacing allows that
    # date to be included.
    pivots = deque([])
    pivot = end
    while pivot >= start:
        pivots.appendleft(pivot)
        pivot = end - delta * len(pivots)

    pivots = np.array([numeric_date(pivot) for pivot in pivots])

    return np.around(pivots, 4)


def make_pivots(pivots, tps):
    '''
    if pivots is a scalar, make a grid of pivot points covering the entire range

    Parameters
    ----------
    pivots : int or iterable
        either number of pivots (a scalar) or the actual pivots
        (will be cast to array and returned)
    tps : numpy.ndarray
        observation time points. Will generate pivots spanning min/max

    Returns
    -------
    pivots : numpy.ndarray
        array of pivot values
    '''
    if isinstance(pivots, int):
        dt = np.max(tps)-np.min(tps)
        return np.linspace(np.min(tps)-0.01*dt, np.max(tps)+0.01*dt, pivots)
    else:
        return np.array(pivots)


def count_observations(pivots, tps):
    pivots = make_pivots(pivots, tps)
    dt = pivots[1]-pivots[0]
    counts, _ = np.histogram(tps, bins=[pivots[0]-0.5*dt] + list(pivots+0.5*dt))
    return counts

def running_average(obs, ws):
    '''
    calculates a running average
    obs     --  observations
    ws      --  window size (number of points to average)

    Parameters
    ----------
    obs : list or numpy.ndarray(bool)
        observations
    ws : int
        window size as measured in number of consecutive points

    Returns
    -------
    numpy.ndarray(float)
        running average of the boolean observations
    '''
    ws=int(ws)
    try:
        tmp_vals = np.convolve(np.ones(ws, dtype=float)/ws, obs, mode='same')
        # fix the edges. using mode='same' assumes zeros outside the range
        if ws%2==0:
            tmp_vals[:ws//2]*=float(ws)/np.arange(ws//2,ws)
            if ws//2>1:
                tmp_vals[-ws//2+1:]*=float(ws)/np.arange(ws-1,ws//2,-1.0)
        else:
            tmp_vals[:ws//2]*=float(ws)/np.arange(ws//2+1,ws)
            tmp_vals[-ws//2:]*=float(ws)/np.arange(ws,ws//2,-1.0)
    except:
        tmp_vals = 0.5*np.ones_like(obs, dtype=float)
    return tmp_vals


def fix_freq(freq, pc):
    '''
    restricts frequencies to the interval [pc, 1-pc]
    removes np.nan values and avoids taking logarithms of 0 or divisions by 0

    Parameters
    ----------
    freq : numpy.ndarray
        frequency trajectory to be thresholded
    pc : float
        threshold value

    Returns
    -------
    numpy.ndarray
        thresholded frequency trajectory
    '''
    freq[np.isnan(freq)]=pc
    return np.minimum(1-pc, np.maximum(pc,freq))


def logit_transform(freq, pc):
    return np.log(np.maximum(freq, pc)/np.maximum(pc,(1-freq)))


def logit_inv(logit_freq, pc):
    tmp_freq = np.maximum(pc, np.minimum(1.0/pc, np.exp(logit_freq)))
    return tmp_freq/(1.0+tmp_freq)


def pq(p):
    return p*(1-p)


class frequency_estimator(object):
    '''
    estimates a smooth frequency trajectory given a series of time stamped
    0/1 observations. The most likely set of frequencies at specified pivot values
    is determined by numerical minimization. Likelihood consist of a bernoulli sampling
    term as well as a term penalizing rapid frequency shifts. this term is motivated by
    genetic drift, i.e., sampling variation.
    '''

    def __init__(self, tps, obs, pivots, stiffness = 20.0,
                inertia = 0.0,  tol=1e-3, pc=1e-4, ws=100,
                method='powell', **kwargs):
        """basic frequency estimator class

        Parameters
        ----------
        tps : list or numpy.ndarray(float)
            array with numerical dates
        obs : list or numpy.ndarray(bool)
            array with boolean observations
        pivots : int or numpy.ndarray(float)
            either integer specifying the number of pivot values,
            or list of explicity pivots
        stiffness : float, optional
            parameter determining how much rapid changes in frequency
            are penalized
        inertia : float, optional
            parameter specifying the prior for the frequency derivitative.
            if intertia=1, the prior is that the slope doesn't change.
            if intertia=0, the prior is that the slope is zero
        tol : float, optional
            optimization tolerance
        pc : float, optional
            pseudo-count/minimal frequency
        ws : int, optional
            window size used to estimate the initial guess for the
            frequency trajectory
        method : str, optional
            optimization method passed down to scipy.minimize
        **kwargs
            Description
        """
        tmp_obs = np.array(sorted(zip(tps, obs), key=lambda x:x[0]))
        self.tps = tmp_obs[:,0]
        self.obs = np.array(tmp_obs[:,1], dtype=bool)
        self.stiffness = stiffness

        self.inertia = inertia
        self.interpolation_type = 'linear'
        self.tol = tol
        self.reg = 1e-6
        self.pc = pc
        self.ws = ws
        self.verbose = 0
        self.method = method

        self.pivots = make_pivots(pivots, self.tps)

        good_tps = (self.tps>=self.pivots[0])&(self.tps<self.pivots[-1])
        self.tps = self.tps[good_tps]
        self.obs = self.obs[good_tps]


    def initial_guess(self, pc=0.01):
        # generate a useful initial guess from a running average of the counts
        if self.ws<len(self.obs):
            tmp_vals = running_average(self.obs, self.ws)
        else:
            tmp_vals = running_average(self.obs, len(self.obs))

        tmp_interpolator = interp1d(self.tps, tmp_vals, bounds_error=False, fill_value = -1)
        pivot_freq = tmp_interpolator(self.pivots)
        pivot_freq[self.pivots<=tmp_interpolator.x[0]] = tmp_vals[0]
        pivot_freq[self.pivots>=tmp_interpolator.x[-1]] = tmp_vals[-1]
        pivot_freq =fix_freq(pivot_freq, pc)
        return pivot_freq

    def stiffLH(self):
        freq = self.pivot_freq
        dfreq = np.diff(freq)
        dfreq_penalty = dfreq[1:] - self.inertia*dfreq[:-1]
        pq_dt_freq = self.dt * pq(freq[:-1])

        # return wright fisher diffusion likelihood for frequency change.
        return -0.25*self.stiffness*\
                    (np.sum(dfreq_penalty**2/pq_dt_freq[1:])+dfreq[0]**2/pq_dt_freq[0])


    def learn(self, initial_guess=None):
        # print('optimizing with', len(self.pivots), 'pivots')
        self.dt = np.diff(self.pivots)
        def logLH(x):
            # x is logit frequency x = log(p/(1-p)) -> p = e^x/(1+e^x); 1-p = 1/(1+e^x)
            self.pivot_freq = logit_inv(x, self.pc)
            try:
                freq = interp1d(self.pivots, x, kind=self.interpolation_type,bounds_error = False, assume_sorted=True) # Throws a bug with some numpy installations, isn't necessary other than for speed.
            except:
                freq = interp1d(self.pivots, x, kind=self.interpolation_type,bounds_error = False)
            estfreq = freq(self.tps)
            # log(p^obs (1-p)^(1-obs)) = log((p/(1-p))^obs (1-p)) = obs*x - log(1+e^x)
            bernoulli_LH = np.sum(estfreq[self.obs]) - np.sum(np.log(1+np.exp(estfreq)))

            stiffness_LH = self.stiffLH()
            LH = stiffness_LH + bernoulli_LH
            return -LH + np.sum((np.abs(x)>log_thres)*x**2)

        from scipy.optimize import minimize
        if initial_guess is None:
            initial_freq = self.initial_guess(pc=0.01)
        else:
            initial_freq = initial_guess(self.pivots)

        self.frequency_estimate = interp1d(self.pivots, initial_freq, kind=self.interpolation_type, bounds_error=False)

        self.pivot_freq = self.frequency_estimate(self.pivots)

        # determine the optimal pivot freqquencies
        self.sol = minimize(logLH, logit_transform(self.pivot_freq, pc=self.pc), method=self.method)
        if self.sol['success']:
            self.pivot_freq = logit_inv(self.sol['x'], self.pc)
            #self.covariance = self.sol['hess_inv'].todense()
            #self.conf95 = 1.96*np.sqrt(np.diag(self.covariance))
        else:
            if debug: import ipdb; ipdb.set_trace()
            if self.method != "powell":
                print("Optimization failed, trying with powell")
                self.sol = minimize(logLH, logit_transform(self.pivot_freq,self.pc), method='powell')

            self.pivot_freq = logit_inv(self.sol['x'], self.pc)
        self.pivot_freq = fix_freq(self.pivot_freq, self.pc)

        # instantiate an interpolation object based on the optimal frequency pivots
        self.frequency_estimate = interp1d(self.pivots, self.pivot_freq, kind=self.interpolation_type, bounds_error=False)

        if self.verbose: print ("neg logLH using",len(self.pivots),"pivots:", self.logLH(self.pivot_freq))


class freq_est_clipped(object):
    """
    simple wrapper for the frequency estimator that attempts
    to estimate a frequency trajectory on a sensible range of
    pivots to avoid frequency optimization at points without data

    Attributes
    ----------
    dtps
        Description
    fe
        Description
    good_pivots
        Description
    good_tps
        Description
    obs
        Description
    pivot_freq
        Description
    pivot_lower_cutoff
        Description
    pivot_upper_cutoff
        Description
    pivots
        Description
    tps
        Description
    valid : bool
        Description
    """
    def __init__(self, tps, obs, pivots, dtps=None, **kwargs):
        super(freq_est_clipped, self).__init__()
        tmp_obs = np.array(sorted(zip(tps, obs), key=lambda x:x[0]))
        self.tps = tmp_obs[:,0]
        self.obs = np.array(tmp_obs[:,1], dtype=bool)
        self.pivots = pivots
        pivot_dt = self.pivots[1]-self.pivots[0]
        if dtps==None:
            self.dtps = 6.0*pivot_dt
        else:
            self.dtps = max(dtps, pivot_dt)

        cum_obs = np.diff(self.obs).cumsum()
        first_obs = max(pivots[0], self.tps[cum_obs.searchsorted(cum_obs[0]+1)])
        last_obs = min(pivots[-1], max(first_obs,self.tps[min(len(self.tps)-1, 20+cum_obs.searchsorted(cum_obs[-1]))]))
        tps_lower_cutoff = first_obs - self.dtps
        tps_upper_cutoff = last_obs + self.dtps


        self.good_tps = (self.tps>=tps_lower_cutoff)&(self.tps<tps_upper_cutoff)&(self.tps<self.pivots[-1])&(self.tps>=self.pivots[0])
        self.valid = True
        if self.good_tps.sum()<3:
            print("too few valid time points:", self.good_tps.sum())
            self.valid=False
            return None
        if self.good_tps.sum()<7:
            from scipy.ndimage import binary_dilation
            self.good_tps = binary_dilation(self.good_tps, iterations=5)
        reduced_obs = self.obs[self.good_tps]
        reduced_tps = self.tps[self.good_tps]
        self.pivot_lower_cutoff = min(reduced_tps[0], tps_lower_cutoff)-pivot_dt
        self.pivot_upper_cutoff = max(reduced_tps[-1], tps_upper_cutoff)+pivot_dt

        self.good_pivots = (self.pivots>=self.pivot_lower_cutoff)\
                             &(self.pivots<self.pivot_upper_cutoff)
        if self.good_pivots.sum()<2:
            from scipy.ndimage import binary_dilation
            self.good_pivots = binary_dilation(self.good_pivots, iterations=2)

        self.fe = frequency_estimator(reduced_tps, reduced_obs,
                                  self.pivots[self.good_pivots], **kwargs)


    def learn(self):
        self.fe.learn()

        self.pivot_freq = np.zeros_like(self.pivots)
        self.pivot_freq[self.good_pivots] = self.fe.pivot_freq

        # set pivots outside of the window used for estimation to strictly zero or one
        # we previously had simply filled all data points with closest pivot that was
        # estimated, but this potentially causes propagation of numerical issues when
        # estimating many nested clades.
        if self.fe.pivot_freq[0]<0.5:
            self.pivot_freq[self.pivots<self.pivot_lower_cutoff] = 0.0
        else:
            self.pivot_freq[self.pivots<self.pivot_lower_cutoff] = 1.0

        if self.fe.pivot_freq[-1]<0.5:
            self.pivot_freq[self.pivots>=self.pivot_upper_cutoff] = 0.0
        else:
            self.pivot_freq[self.pivots>=self.pivot_upper_cutoff] = 1.0


class nested_frequencies(object):
    """
    estimates frequencies of mutually exclusive events such as mutations
    at a particular genomic position or subclades in a tree
    """
    def __init__(self, tps, obs, pivots,  **kwargs):
        """Summary

        Parameters
        ----------
        tps : numpy.ndarray
            array of numerical dates
        obs : numpy.ndarray(bool)
            array of true/false observations
        pivots : numpy.ndarray
            pivot values
        **kwargs
            Description
        """
        super(nested_frequencies, self).__init__()
        self.tps = tps
        self.obs = obs
        self.pivots = pivots
        self.kwargs = kwargs

    def calc_freqs(self):
        sorted_obs = sorted(self.obs.items(), key=lambda x:x[1].sum(), reverse=True)
        self.remaining_freq = np.ones_like(self.pivots)

        self.frequencies = {}
        valid_tps = np.ones_like(self.tps, dtype=bool)
        for mut, obs in sorted_obs[:-1]:
            # print(mut,'...')
            fe = freq_est_clipped(self.tps[valid_tps], obs[valid_tps], self.pivots, **self.kwargs)
            if fe.valid==False:
                self.frequencies[mut] = np.zeros_like(self.remaining_freq)
                break

            fe.learn()
            self.frequencies[mut] = self.remaining_freq * fe.pivot_freq
            self.remaining_freq *= (1.0-fe.pivot_freq)
            valid_tps = valid_tps&(~obs)

        self.frequencies[sorted_obs[-1][0]] = self.remaining_freq
        return self.frequencies


class tree_frequencies(object):
    '''
    class that estimates frequencies for nodes in the tree. each internal node is assumed
    to be named with an attribute clade, of root doesn't have such an attribute, clades
    will be numbered in preorder. Each node is assumed to have an attribute `attr` with a
    key "num_date".
    '''
    def __init__(self, tree, pivots, node_filter=None, min_clades=10, verbose=0, pc=1e-4, **kwargs):
        '''
        set up the internal tree, the pivots and cutoffs

        Parameters
        ----------
        tree : Bio.Phylo.BaseTree.Tree
            Biopython tree
        pivots : int or array
            number or list of pivots
        node_filter : callable, optional
            function that evaluates to true/false to filter nodes
        min_clades : int, optional
            minimal size of clades to estimate proper frequencies
        verbose : int, optional
            verbosity
        pc : float, optional
            pseudo-counts/lower cutoff to drive estimates away from 0/1
        **kwargs
            Description
        '''
        self.tree = tree
        self.min_clades = min_clades
        self.pivots = pivots
        self.kwargs = kwargs
        self.verbose = verbose
        self.pc = pc
        if node_filter is None:
            self.node_filter = lambda x:True
        else:
            self.node_filter = node_filter
        self.prepare()


    def prepare(self):
        # name nodes if they aren't already named
        if not hasattr(self.tree.root, 'clade'):
            for ni,node in enumerate(self.tree.find_clades(order='preorder')):
                node.clade = ni

        # extract all observations
        tps = []
        leaf_count = 0
        for node in self.tree.find_clades(order='postorder'):
            if node.is_terminal():
                if self.node_filter(node):
                    tps.append(node.attr["num_date"])
                    node.leafs = np.array([leaf_count], dtype=int)
                    leaf_count+=1
                else:
                    node.leafs = np.array([], dtype=int)
            else:
                node.leafs = np.concatenate([c.leafs for c in node.clades])
        self.tps = np.array(tps)

        if np.isscalar(self.pivots):
            self.frequencies = {self.tree.root.clade:np.ones(self.pivots)}
        else:
            self.frequencies = {self.tree.root.clade:np.ones_like(self.pivots)}

        self.counts = count_observations(self.pivots, self.tps)


    def estimate_clade_frequencies(self):
        for node in self.tree.get_nonterminals(order='preorder'):
            if self.verbose: print("Estimating frequencies of children of node ",node.clade)
            node_tps = self.tps[node.leafs]
            obs_to_estimate = {}
            small_clades = []
            if len(node.clades)==1:
                self.frequencies[node.clades[0].clade] = self.frequencies[node.clade]
                continue

            for c in node.clades:
                if len(c.leafs)>self.min_clades:
                    obs_to_estimate[c.clade] = np.in1d(node.leafs, c.leafs)
                else:
                    # Only include internal nodes or tips that pass the node
                    # filter as small clades.
                    if not c.is_terminal() or self.node_filter(c):
                        small_clades.append(c)
            if len(obs_to_estimate):
                if len(small_clades):
                    remainder = {}
                    for c in small_clades:
                        remainder[c.clade] = np.in1d(node.leafs, c.leafs)
                    if len(small_clades)==1:
                        obs_to_estimate.update(remainder)
                    else:
                        obs_to_estimate['other'] = np.any(list(remainder.values()), axis=0)

                ne = nested_frequencies(node_tps, obs_to_estimate, self.pivots, pc=self.pc, **self.kwargs)
                freq_est = ne.calc_freqs()
                for clade, tmp_freq in freq_est.items():
                    if clade != "other":
                        self.frequencies[clade] = self.frequencies[node.clade] * tmp_freq

                if len(small_clades) > 1:
                    total_leaves_in_small_clades = 0
                    for clade in small_clades:
                        total_leaves_in_small_clades += len(clade.leafs)

                    for clade in small_clades:
                        if total_leaves_in_small_clades > 0:
                            frac = len(clade.leafs) / total_leaves_in_small_clades
                        else:
                            frac = 0.0

                        self.frequencies[clade.clade] = frac * freq_est["other"] * self.frequencies[node.clade]
            else:
                for clade in small_clades:
                    if len(node.leafs):
                        frac = 1.0*len(clade.leafs)/len(node.leafs)
                    else:
                        frac = 0.0
                    self.frequencies[clade.clade] = frac*self.frequencies[node.clade]

        # Assign zero frequencies to tips that did not pass the node_filter.
        for tip in self.tree.get_terminals():
            if not self.node_filter(tip):
                self.frequencies[tip.clade] = np.zeros(len(self.pivots))


    def calc_confidence(self):
        '''
        for each frequency trajectory, calculate the bernouilli sampling error
        -- in reality we should look at the inverse hessian, but this is a
        useful approximation in most cases

        Returns
        -------
        dict
            dictionary with estimated confidence intervals
        '''
        self.confidence = {}
        for key, freq in self.frequencies.items():
            # add a pseudo count 1/(n+1) and normalize to n+1
            self.confidence[key] = np.sqrt((1.0/(1+self.counts)+freq*(1-freq))/(1.0+self.counts))
        return self.confidence


class alignment_frequencies(object):
    '''
    calculates frequencies of mutations in an alignment. uses nested frequencies
    for mutations at a particular site such that mutations are forced to add
    up to one.

    '''
    def __init__(self, aln, tps, pivots, **kwargs):
        """Create an instance of the alignment frequency estimate

        Parameters
        ----------
        aln : Bio.Align.MultipleSeqAlignment
            alignment
        tps : np.ndarray(float)
            Array of numerical dates, one for each sequence in the
            alignment in the SAME ORDER!
        pivots : np.ndarray(float)
            pivot values for which frequencies are estimated
        **kwargs
            Description
        """
        self.aln = np.array(aln)
        self.tps = np.array(tps)
        self.kwargs = kwargs
        self.pivots = pivots
        self.counts = count_observations(self.pivots, self.tps)


    def estimate_genotype_frequency(self, aln, gt, **kwargs):
        '''
        slice an alignment at possibly multiple positions and calculate the
        frequency trajectory of this multi-locus genotype

        Parameters
        ----------
        gt : list
            a list of (position, state) tuples specifying the genotype
            whose frequency is to be estimated

        Returns
        -------
        numpy.ndarray
            frequency trajectory
        '''
        match = []
        for pos, state in gt:
            match.append(aln[:,pos]==state)
        obs = np.array(match).all(axis=0)

        fe = frequency_estimator(zip(self.tps, obs),self.pivots, **kwargs)
        fe.learn()
        return fe.frequency_estimate


    def mutation_frequencies(self, min_freq=0.01, include_set=None, ignore_char=''):
        '''
        estimate frequencies of single site mutations for each alignment column.
        This function populates a dictionary class.frequencies with the frequency
        trajectories.

        Parameters
        ----------
        min_freq : float, optional
            minimal all-time frequency for an aligment column to be considered
        include_set : list or set, optional
            set of alignment column that will be used regardless of variation
        ignore_char : str, optional
            ignore this character in an alignment column (missing data)
        '''
        if include_set is None:
            include_set=[]
        alphabet = np.unique(self.aln)
        # af: rows correspond to letters of the alphabet, columns positions in the alignemtn, values: frequencies
        af = np.zeros((len(alphabet), self.aln.shape[1]))
        for ni, n in enumerate(alphabet):
            af[ni] = (self.aln==n).mean(axis=0)

        # minor_freqs: frequencies of all alleles not the major one! ignore_char is excluded
        if ignore_char:
            minor_freqs = af[alphabet!=ignore_char].sum(axis=0) - af[alphabet!=ignore_char].max(axis=0)
        else:
            minor_freqs = 1.0 - af.max(axis=0)
        self.frequencies = {}
        for pos in sorted(set.union(set(np.where(minor_freqs>min_freq)[0]), include_set)):
            # indices determine the indicies of mutations in descending order by frequency
            nis = np.argsort(af[:,pos])[::-1]
            nis = nis[af[nis,pos]>0]
            column = self.aln[:,pos]
            if ignore_char: # subset sequences, time points and alphabet to non-gapped and non X
                good_seq = (column!=ignore_char)&(column!='X') #array of bools, length = # seqs in aln
                tps=self.tps[good_seq]
                column = column[good_seq]
                nis = nis[(alphabet[nis]!=ignore_char)&(alphabet[nis]!='X')]
            else:
                tps = self.tps

            if len(nis)==0:
                # no characters with any frequencies (column of ignore chars)
                continue

            muts = alphabet[nis]
            obs = {}
            for ni, mut in zip(nis, muts):
                if (af[ni,pos]>min_freq and af[ni,pos]<1-min_freq) or ni==0 or (pos in include_set):
                    obs[(pos, mut)] = column==mut
                else:
                    break

            # if multiple states are below frequency threshold, glob them together as 'other'
            if len(obs)==0:
                obs[(pos, muts[0])] = column==muts[0]
            elif len(obs)!=len(nis):
                tmp = ~np.any(list(obs.values()), axis=0) # pull out sequences not yet assigned

                if any(tmp):
                    if len(obs)==len(nis)-1: # if only category left, assign it
                        obs[(pos, muts[-1])] = tmp
                    else: #other wise, call that mutation 'other'
                        obs[(pos, 'other')] = tmp

            print("Estimating frequencies of position: {}".format(pos), end="\r")
            sys.stdout.flush()
            # print("Variants found at frequency:", [(k,o.mean()) for k,o in obs.items()])

            # calculate frequencies, which will be added to the frequencies dict with (pos, mut) as key
            ne = nested_frequencies(tps, obs, self.pivots, **self.kwargs)
            self.frequencies.update(ne.calc_freqs())


    def calc_confidence(self):
        """calculate a crude binomial sampling confidence interval
        of the frequency estimate. This ignores autocorrelation of the
        trajectory and returns one standard deviation (68%).

        Returns
        -------
        dict
            dictionary of standard deviations for each estimated trajectory
        """
        self.confidence = {}
        for key, freq in self.frequencies.items():
            # add a pseudo count 1/(n+1) and normalize to n+1
            self.confidence[key] = np.sqrt((1.0/(1+self.counts)+freq*(1-freq))/(1.0+self.counts))

        return self.confidence


def test_simple_estimator():
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        plot=False
        print("Plotting requires a working matplotlib installation.")

    tps = np.sort(100 * np.random.uniform(size=500))
    freq_traj = [0.1]
    stiffness=100
    s=-0.02
    for dt in np.diff(tps):
        freq_traj.append(freq_traj[-1]*np.exp(-s*dt)+np.sqrt(2*np.max(0,freq_traj[-1]*(1-freq_traj[-1]))*dt/stiffness)*np.random.normal())
    obs = np.random.uniform(size=tps.shape)<np.array(freq_traj)

    fe = frequency_estimator(tps, obs, pivots=20, stiffness=stiffness, inertia=0.7)
    fe.learn()
    freq = fe.frequency_estimate(fe.tps)

    if plot:
        plt.figure()
        plt.plot(tps, freq_traj, 'o', label = 'actual frequency')
        plt.plot(fe.tps, freq, '-', label='interpolation')
        plt.plot(tps, (2*obs-1)*0.05, 'o')
        plt.plot(tps[obs], 0.05*np.ones(np.sum(obs)), 'o', c='r', label = 'observations')
        plt.plot(tps[~obs], -0.05*np.ones(np.sum(1-obs)), 'o', c='g')
        plt.plot(tps, np.zeros_like(tps), 'k')
        ws=20
        r_avg = running_average(obs, ws)
        plt.plot(fe.tps[ws/2:-ws/2+1], np.convolve(np.ones(ws, dtype=float)/ws, obs, mode='valid'), 'r', label = 'running avg')
        plt.plot(fe.tps, r_avg, 'k', label = 'running avg')
        plt.legend(loc=2)
    return fe

def test_nested_estimator():
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        plot=False
        print("Plotting requires a working matplotlib installation.")

    tps = np.sort(100 * np.random.uniform(size=2000))
    freq_traj = [0.1]
    stiffness=1000
    s=-0.03
    for dt in np.diff(tps):
        freq_traj.append(freq_traj[-1]*np.exp(-s*dt)+np.sqrt(2*np.max(0,freq_traj[-1]*(1-freq_traj[-1]))*dt/stiffness)*np.random.normal())
    freq_traj = np.array(freq_traj)
    tmp = np.random.uniform(size=tps.shape)
    obs={}
    obs['A'] = tmp<freq_traj*0.5
    obs['B'] = (tmp>=freq_traj*0.5)&(tmp<freq_traj*0.7)
    obs['C'] = tmp>=freq_traj*0.7

    fe = nested_frequencies(tps, obs, pivots=20, stiffness=stiffness, inertia=0.7)
    nested_freq = fe.calc_freqs()

    if plot:
        plt.figure()
        for k,o in obs.items():
            plt.plot(tps, o, 'o')
            plt.plot(fe.pivots,fe.frequencies[k], 'o')

    return nested_freq

def float_to_datestring(numdate):
    """convert a numeric decimal date to a python datetime object
    Note that this only works for AD dates since the range of datetime objects
    is restricted to year>1.
    Copied from treetime.utils
    Parameters
    ----------
    numdate : float
        numeric date as in 2018.23
    Returns
    -------
    datetime.datetime
        datetime object
    """
    from calendar import isleap
    days_in_year = 366 if isleap(int(numdate)) else 365
    # add a small number of the time elapsed in a year to avoid
    # unexpected behavior for values 1/365, 2/365, etc
    days_elapsed = int(((numdate%1)+1e-10)*days_in_year)
    date = datetime.datetime(int(numdate),1,1) + datetime.timedelta(days=days_elapsed)

    return "%s-%02d-%02d" % (date.year, date.month, date.day)

def timestamp_to_float(time):
    """Convert a pandas timestamp to a floating point date.
    This is not entirely accurate as it doesn't account for months with different
    numbers of days, but should be close enough to be accurate for weekly pivots.

    Examples
    --------
    >>> import datetime
    >>> time = datetime.date(2010, 10, 1)
    >>> timestamp_to_float(time)
    2010.75
    >>> time = datetime.date(2011, 4, 1)
    >>> timestamp_to_float(time)
    2011.25
    >>> timestamp_to_float(datetime.date(2011, 1, 1))
    2011.0
    >>> timestamp_to_float(datetime.date(2011, 12, 1)) == (2011.0 + 11.0 / 12)
    True
    """
    return time.year + ((time.month - 1) / 12.0) + ((time.day - 1) / 365.25)


class KdeFrequencies(object):
    """Methods to estimate clade frequencies for phylogenetic trees by creating
    normal distributions from timestamped tips in the tree and building a kernel
    density estimate across discrete time points from these tip observations for
    each clade in the tree.
    """
    def __init__(self, sigma_narrow=1 / 12.0, sigma_wide=3 / 12.0, proportion_wide=0.2,
                 pivot_frequency=1, start_date=None, end_date=None,
                 pivot_interval_units="months", weights=None, weights_attribute=None,
                 node_filters=None, max_date=None, include_internal_nodes=False,
                 censored=False):
        """Define parameters for KDE-based frequency estimation.

        Args:
            sigma_narrow (float): Bandwidth for first of two Gaussians composing the KDEs
            sigma_wide (float): Bandwidth for second of two Gaussians composing the KDEs
            proportion_wide (float): Proportion of the second Gaussian to include in each KDE
            pivot_frequency (int): Number of months between pivots
            start_date (float): start of the pivots interval
            end_date (float): end of the pivots interval
            pivot_interval_units (str): Whether pivot intervals are measured in "months" or "weeks"
            weights (dict): Numerical weights indexed by attribute values and applied to individual tips
            weights_attribute (str): Attribute annotated on tips of a tree to use for weighting
            node_filters (dict): Mapping of node attribute names (keys) to a list of valid values to keep
            max_date (float): Maximum year beyond which tips are excluded from frequency estimation and are assigned
                              frequencies of zero
            include_internal_nodes (bool): Whether internal (non-tip) nodes should have their frequencies estimated
            censored (bool): Whether future observations should be censored at each pivot

        Returns:
            KdeFrequencies
        """
        self.sigma_narrow = sigma_narrow
        self.sigma_wide = sigma_wide
        self.proportion_wide = proportion_wide
        self.pivot_frequency = pivot_frequency
        self.start_date = start_date
        self.end_date = end_date
        self.pivot_interval_units = pivot_interval_units
        self.weights = weights
        self.weights_attribute = weights_attribute
        self.node_filters = node_filters
        self.max_date = max_date
        self.include_internal_nodes = include_internal_nodes
        self.censored = censored

    def get_params(self):
        """
        Returns the parameters used to define the current instance.

        Returns:
            dict: parameters that define the current instance and that can be used to create a new instance
        """
        return {
            "sigma_narrow": self.sigma_narrow,
            "sigma_wide": self.sigma_wide,
            "proportion_wide": self.proportion_wide,
            "pivot_frequency": self.pivot_frequency,
            "start_date": self.start_date,
            "end_date": self.end_date,
            "pivot_interval_units": self.pivot_interval_units,
            "weights": self.weights,
            "weights_attribute": self.weights_attribute,
            "max_date": self.max_date,
            "include_internal_nodes": self.include_internal_nodes,
            "node_filters": self.node_filters,
            "censored": self.censored
        }

    @classmethod
    def from_json(cls, json_dict):
        """Returns an instance populated with parameters and data from the given JSON dictionary.
        """
        params = json_dict["params"]
        instance = cls(**params)

        if "data" in json_dict:
            instance.pivots = np.array(json_dict["data"]["pivots"])
            frequencies = json_dict["data"]["frequencies"]

            instance.frequencies = {}
            for clade in frequencies:
                instance.frequencies[clade] = np.array(frequencies[clade])

        return instance

    def to_json(self):
        """Returns a dictionary for the current instance that can be serialized in a JSON file.
        """
        frequencies_json = {
            "params": self.get_params()
        }

        # If frequencies have been estimated, export them along with the pivots as data.
        if hasattr(self, "frequencies"):
            frequencies = {}
            for clade in self.frequencies:
                # numpy arrays are not supported by JSON and need to be converted to lists.
                frequencies[clade] = self.frequencies[clade].tolist()

            frequencies_json["data"] = {
                "pivots": self.pivots.tolist(),
                "frequencies": frequencies
            }

        return frequencies_json

    @classmethod
    def get_density_for_observation(cls, mu, pivots, sigma_narrow=1/12.0, sigma_wide=3/12.0, proportion_wide=0.2, **kwargs):
        """Build a normal distribution centered across the given floating point date,
        mu, with a standard deviation based on the given sigma value and return
        the probability mass at each pivot. These mass values per pivot will form the
        input for a kernel density estimate across multiple observations.
        """
        return ((1-proportion_wide) * norm.pdf(pivots, loc=mu, scale=sigma_narrow) +
                proportion_wide * norm.pdf(pivots, loc=mu, scale=sigma_wide))

    @classmethod
    def get_densities_for_observations(cls, observations, pivots, max_date=None, **kwargs):
        """Create a matrix of densities for one or more observations across the given
        pivots with one row per observation and one column per pivot.

        Observations can be optionally filtered by a maximum date such that all
        densities are estimated to be zero after that date.
        """
        density_matrix = np.zeros((len(observations), len(pivots)))
        for i, obs in enumerate(observations):
            # If the given observation occurred outside the given pivots or
            # after the given maximum date, do not try to estimate its
            # frequency.
            if (obs < pivots[0] or obs > pivots[-1]) or (max_date is not None and obs > max_date):
                density = np.zeros_like(pivots)
            else:
                density = cls.get_density_for_observation(obs, pivots, **kwargs)

            density_matrix[i] = density

        return density_matrix

    @classmethod
    def normalize_to_frequencies(cls, density_matrix, normalize_to=1.0):
        """Normalize the values of a given density matrix to 1 across all columns
        (time points) with non-zero sums. This converts kernal PDF mass into a
        frequency estimate.
        """
        normalized_freq_matrix = density_matrix.copy()

        # Find columns that can be meaningfully normalized.
        nonzero_columns = np.where(density_matrix.sum(axis=0) > 0)[0]

        # Normalize by column.
        normalized_freq_matrix[:, nonzero_columns] = (normalize_to * density_matrix[:, nonzero_columns] /
                                                      density_matrix[:, nonzero_columns].sum(axis=0))

        return normalized_freq_matrix

    @classmethod
    def estimate_frequencies(cls, tip_dates, pivots, normalize_to=1.0, max_date=None, **kwargs):
        """Estimate frequencies of the given observations across the given pivots.
        """
        # Calculate base frequencies from observations.
        if kwargs.get("censored"):
            # Calculate censored frequencies at each pivot. If a maximum date
            # has already been requested, only calculate censored frequencies up
            # to that point.
            if max_date is None:
                max_date = pivots[-1]

            pivots_to_censor = [pivot for pivot in pivots if pivot <= max_date]
            density_matrix = None
            for i, pivot in enumerate(pivots_to_censor):
                # Censor observations from the future using the smallest value
                # of either the current pivot or the requested max date.
                censored_matrix = cls.get_densities_for_observations(
                    tip_dates,
                    pivots,
                    max_date=min(pivot, max_date),
                    **kwargs
                )

                if density_matrix is None:
                    density_matrix = censored_matrix
                else:
                    density_matrix[:, i] = censored_matrix[:, i]
        else:
            density_matrix = cls.get_densities_for_observations(tip_dates, pivots, max_date=max_date, **kwargs)

        # Normalize frequencies to sum to 1.
        normalized_freq_matrix = cls.normalize_to_frequencies(density_matrix, normalize_to=normalize_to)

        return normalized_freq_matrix


class TreeKdeFrequencies(KdeFrequencies):
    """KDE frequency estimator for phylogenetic trees

    Estimates frequencies for samples provided as annotated tips on a given
    tree and optionally reports clade-level frequencies based on the sum of each
    clade's respective tips.
    """
    def tip_passes_filters(self, tip):
        """Returns a boolean indicating whether a given tip passes the node filters
        defined for the current instance.

        If no filters are defined, returns True.

        Args:
            tip (Bio.Phylo.BaseTree.Tree): tip from a Bio.Phylo tree annotated with attributes in `tip.attr`

        Returns:
            bool: whether the given tip passes the defined filters or not
        """
        return (self.node_filters is None or
                all([tip.attr[key] in values for key, values in self.node_filters.items()]))

    def estimate_tip_frequencies_to_proportion(self, tips, proportion):
        """Estimate frequencies for a given set of tips up to a given proportion of total frequencies.

        Args:
            tips (list): a list of Bio.Phylo terminals annotated with attributes in `tip.attr`
            proportion (float): the proportion of the total frequency that the given tips should represent

        Returns:
            dict: frequencies of given tips indexed by tip name
        """
        clade_frequencies = {}

        # If none of the tips pass the given node filters, do not try to
        # normalize tip frequencies.
        if len(tips) == 0:
            return clade_frequencies

        tips = np.array(sorted(tips, key=lambda row: row[1]))
        clades = tips[:, 0]
        tip_dates = tips[:, 1].astype(float)

        # Map clade ids to their corresponding frequency matrix row index.
        clade_to_index = {clades[i]: i for i in range(len(clades))}

        # Calculate tip frequencies and normalize.
        normalized_freq_matrix = self.estimate_frequencies(
            tip_dates,
            self.pivots,
            normalize_to=proportion,
            sigma_narrow=self.sigma_narrow,
            sigma_wide=self.sigma_wide,
            proportion_wide=self.proportion_wide,
            max_date=self.max_date,
            censored=self.censored
        )

        for clade in clades:
            clade_frequencies[clade] = normalized_freq_matrix[clade_to_index[clade]]

        return clade_frequencies

    def estimate(self, tree):
        """Estimate frequencies for a given tree using the parameters defined for this instance.

        If weights are defined, frequencies represent a weighted mean across the
        values in attribute defined by `self.weights_attribute`.

        Args:
            tree (Bio.Phylo.BaseTree.Tree): annotated tree whose nodes all have an `attr` attribute with at least  "num_date" key

        Returns:
            dict: node frequencies by clade

        """
        # Calculate pivots for the given tree.
        observations = [tip.attr["num_date"] for tip in tree.get_terminals()]
        self.pivots = get_pivots(
            observations,
            self.pivot_frequency,
            start_date=self.start_date,
            end_date=self.end_date,
            pivot_interval_units=self.pivot_interval_units
        )

        # If weights are defined, calculate frequencies for all tips by their
        # weight attribute values such that the total frequency of all tips with
        # a given value sums to the proportion of tips that have that value.
        # If weights are not defined, estimate frequencies such that they sum to 1.
        if self.weights:
            # Determine which weight attributes are represented in the
            # tree. Drop any unrepresented attributes and renormalize the
            # remaining proportions per attribute to sum to one before
            # estimating frequencies.
            weights_represented = set([
                tip.attr[self.weights_attribute]
                for tip in tree.find_clades(terminal=True)
            ])

            if len(self.weights) != len(weights_represented):
                # Remove unrepresented weights.
                weights_unrepresented = set(self.weights.keys()) - weights_represented
                for weight in weights_unrepresented:
                    del self.weights[weight]

                # Renormalize the remaining weights.
                weight_total = sum(self.weights.values())
                for key, value in self.weights.items():
                    self.weights[key] = value / weight_total

            # Confirm that one or more weights are represented by tips in the
            # tree. If there are no more weights, raise an exception because
            # this likely represents a data error (either in the tree
            # annotations or the weight definitions).
            if len(self.weights) == 0:
                raise TreeKdeFrequenciesError("None of the provided weight attributes were represented by tips in the given tree. Doublecheck weight attribute definitions and their representations in the tree.")

            # Estimate frequencies for all tips within each weight attribute
            # group.
            weight_keys, weight_values = zip(*sorted(self.weights.items()))
            proportions = np.array(weight_values) / np.array(weight_values).sum(axis=0)
            frequencies = {}

            for (weight_key, proportion) in zip(weight_keys, proportions):
                # Find tips with the current weight attribute.
                tips = [(tip.name, tip.attr["num_date"])
                        for tip in tree.get_terminals()
                        if tip.attr[self.weights_attribute] == weight_key and self.tip_passes_filters(tip)]
                frequencies.update(self.estimate_tip_frequencies_to_proportion(tips, proportion))
        else:
            tips = [(tip.name, tip.attr["num_date"])
                    for tip in tree.get_terminals()
                    if self.tip_passes_filters(tip)]
            frequencies = self.estimate_tip_frequencies_to_proportion(tips, proportion=1.0)

        # Assign zero frequencies to any tips that were filtered out of the frequency estimation.
        for tip in tree.get_terminals():
            if not tip.name in frequencies:
                frequencies[tip.name] = np.zeros_like(self.pivots)

        if self.include_internal_nodes:
            for node in tree.find_clades(order="postorder"):
                if not node.is_terminal():
                    frequencies[node.name] = np.array(
                        [frequencies[child.name] for child in node.clades]
                    ).sum(axis=0)

        # Store frequencies in the current instance for simplified exporting.
        self.frequencies = frequencies

        return frequencies


class AlignmentKdeFrequencies(KdeFrequencies):
    """KDE frequency estimator for multiple sequence alignments

    Estimates frequencies for samples provided as sequences in a multiple
    sequence alignment and corresponding observation dates for each sequence.
    """
    def estimate(self, alignment, observations):
        """Estimate frequencies of bases/residues at each site in a given multiple sequence
        alignment based on the observation dates associated with each sequence
        in the alignment.
        """
        # Calculate pivots for the given tree.
        self.pivots = get_pivots(
            observations,
            self.pivot_frequency,
            start_date=self.start_date,
            end_date=self.end_date,
            pivot_interval_units=self.pivot_interval_units
        )

        # Pair alignment sequence indices with observation dates.
        samples = zip(range(len(alignment)), observations)

        # Sort samples by date and extract sequence indices and dates.
        samples = np.array(sorted(samples, key=lambda row: row[1]))
        sample_indices = samples[:, 0]
        sample_dates = samples[:, 1].astype(float)

        # Map sample alignment indices to their corresponding frequency matrix row indices.
        sample_to_index = {sample_indices[i]: i for i in range(len(samples))}

        # Calculate tip frequencies and normalize.
        normalized_freq_matrix = self.estimate_frequencies(
            sample_dates,
            self.pivots,
            sigma_narrow=self.sigma_narrow,
            sigma_wide=self.sigma_wide,
            proportion_wide=self.proportion_wide,
            max_date=self.max_date,
            censored=self.censored
        )

        # Calculate frequencies per site and base/residue from the estimated
        # frequencies per sample.
        frequencies = {}
        for position in range(alignment.get_alignment_length()):
            # Only estimate frequencies at sites with variation.
            if len(set(alignment[:, position])) > 1:
                for sample, base in enumerate(alignment[:, position]):
                    # Estimate frequencies by current position and base/residue
                    # at the position for each sample.
                    key = "%s:%s" % (position + 1, base.upper())
                    if key not in frequencies:
                        frequencies[key] = np.zeros_like(self.pivots)

                    frequencies[key] += normalized_freq_matrix[sample_to_index[sample], :]

        self.frequencies = frequencies

        return frequencies


if __name__=="__main__":
    plot=True
    fe = test_nested_estimator()
    #fe = test_simple_estimator()
