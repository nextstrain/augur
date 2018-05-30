# estimates clade frequencies
from __future__ import division, print_function
from collections import defaultdict
from scipy.interpolate import interp1d
from scipy.stats import norm
import time
import numpy as np
import sys

debug = False
log_thres = 10.0

def make_pivots(pivots, tps):
    '''
    if pivots is a scalar, make a grid of pivot points covering the entire range
    '''
    if np.isscalar(pivots):
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
        import ipdb; ipdb.set_trace()
        tmp_vals = 0.5*np.ones_like(obs, dtype=float)
    return tmp_vals


def fix_freq(freq, pc):
    '''
    restricts frequencies to the interval [pc, 1-pc]
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
                inertia = 0.0,  tol=1e-3, pc=1e-4, ws=100, method='powell', **kwargs):
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

        good_tps = (self.tps>self.pivots[0])&(self.tps<self.pivots[-1])
        self.tps = self.tps[good_tps]
        self.obs = self.obs[good_tps]


    def initial_guess(self, pc=0.01):
        # generate a useful initital case from a running average of the counts
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
            self.dtps = np.max(dtps, pivot_dt)

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
    will be numbered in preorder. Each node is assumed to have an attribute "numdate"
    '''
    def __init__(self, tree, pivots, node_filter=None, min_clades = 20, verbose=0, pc=1e-4, **kwargs):
        '''
        set up the internal tree, the pivots and cutoffs
        node_filter -- a function that can be used to exclude terminals nodes
                       from the estimation. primarily meant for estimating region
                       specific frequencues and training fitness models
        min_clades  -- smallest clades for which frequencies are estimated.
        '''
        self.tree = tree
        self.min_clades = 10 #min_clades
        self.pivots = pivots
        self.kwargs = kwargs
        self.verbose = verbose
        self.pc = pc
        if node_filter is None:
            self.node_filter = lambda x:True
        else:
            self.node_filter = node_filter
        print("filter func", self.node_filter)
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
                    tps.append(node.numdate)
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
                        obs_to_estimate['other'] = np.any(remainder.values(), axis=0)

                ne = nested_frequencies(node_tps, obs_to_estimate, self.pivots, pc=self.pc, **self.kwargs)
                freq_est = ne.calc_freqs()
                for clade, tmp_freq in freq_est.iteritems():
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
        '''
        self.confidence = {}
        for key, freq in self.frequencies.iteritems():
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
        self.aln = np.array(aln)
        self.tps = np.array(tps)
        self.kwargs = kwargs
        self.pivots = pivots
        self.counts = count_observations(self.pivots, self.tps)


    def estimate_genotype_frequency(self, gt):
        '''
        slice an alignment at possibly multiple positions and calculate the
        frequency trajectory of this multi-locus genotype
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
        estimate frequencies of single site mutations for each alignment column
        params
            min_freq:       the minimal minor allele frequency for a column to be included
            include_set:    a set of sites that is to be included regardless of frequencies
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
                tmp = ~np.any(obs.values(), axis=0) # pull out sequences not yet assigned

                if any(tmp):
                    if len(obs)==len(nis)-1: # if only category left, assign it
                        obs[(pos, muts[-1])] = tmp
                    else: #other wise, call that mutation 'other'
                        obs[(pos, 'other')] = tmp

            print("Estimating frequencies of position: {}".format(pos), end="\r")
            sys.stdout.flush()
            # print("Variants found at frequency:", [(k,o.mean()) for k,o in obs.iteritems()])

            # calculate frequencies, which will be added to the frequencies dict with (pos, mut) as key
            ne = nested_frequencies(tps, obs, self.pivots, **self.kwargs)
            self.frequencies.update(ne.calc_freqs())


    def calc_confidence(self):
        self.confidence = {}
        for key, freq in self.frequencies.iteritems():
            # add a pseudo count 1/(n+1) and normalize to n+1
            self.confidence[key] = np.sqrt((1.0/(1+self.counts)+freq*(1-freq))/(1.0+self.counts))

        return self.confidence


def test_simple_estimator():
    import matplotlib.pyplot as plt
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
    import matplotlib.pyplot as plt
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
        for k,o in obs.iteritems():
            plt.plot(tps, o, 'o')
            plt.plot(fe.pivots,fe.frequencies[k], 'o')

    return nested_freq


class KdeFrequencies(object):
    """Methods to estimate clade frequencies for phylogenetic trees by creating
    normal distributions from timestamped tips in the tree and building a kernel
    density estimate across discrete time points from these tip observations for
    each clade in the tree.
    """
    def __init__(self):
        pass

    @classmethod
    def get_density_for_observation(cls, mu, pivots, sigmaNarrow=1/12.0, sigmaWide=3/12.0, proportionWide=0.2):
        """Build a normal distribution centered across the given floating point date,
        mu, with a standard deviation based on the given sigma value and return
        the probability mass at each pivot. These mass values per pivot will form the
        input for a kernel density estimate across multiple observations.
        """
        return (1-proportionWide) * norm.pdf(pivots, loc=mu, scale=sigmaNarrow) + proportionWide * norm.pdf(pivots, loc=mu, scale=sigmaWide)

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
        normalized_freq_matrix[:, nonzero_columns] = normalize_to * density_matrix[:, nonzero_columns] / density_matrix[:, nonzero_columns].sum(axis=0)

        return normalized_freq_matrix

    @classmethod
    def estimate_frequencies(cls, tip_dates, pivots, normalize_to=1.0, **kwargs):
        """Estimate frequencies of the given observations across the given pivots.
        """
        # Calculate base frequencies from observations.
        density_matrix = cls.get_densities_for_observations(tip_dates, pivots, **kwargs)

        # Normalize frequencies to sum to 1.
        normalized_freq_matrix = cls.normalize_to_frequencies(density_matrix, normalize_to=normalize_to)

        return normalized_freq_matrix

    @classmethod
    def estimate_frequencies_for_tree(cls, tree, pivots, **kwargs):
        """Estimate global frequencies for all nodes in a tree across the given pivots.
        """

        clade_frequencies = defaultdict(dict)

        # Find tips within region.
        tips = [(tip.clade, tip.attr["num_date"]) for tip in tree.get_terminals()]
        tips = np.array(sorted(tips, key=lambda row: row[1]))
        clades = tips[:, 0].astype(int)
        tip_dates = tips[:, 1].astype(float)

        # Map clade ids to their corresponding frequency matrix row index.
        clade_to_index = {clades[i]: i for i in range(len(clades))}

        # Calculate tip frequencies and normalize.
        normalized_freq_matrix = cls.estimate_frequencies(tip_dates, pivots, normalize_to=1.0, **kwargs)

        for clade in clades:
            clade_frequencies["global"][clade] = normalized_freq_matrix[clade_to_index[clade]]

        for node in tree.find_clades(order="postorder"):
            if not node.is_terminal():
                clade_frequencies["global"][node.clade] = np.array([clade_frequencies["global"][child.clade]
                                                                for child in node.clades]).sum(axis=0)

        return clade_frequencies


    @classmethod
    def estimate_region_weighted_frequencies_for_tree(cls, tree, pivots, regions, weights, **kwargs):
        """Estimate global and regional frequencies for all nodes in a tree across the
        given pivots. Global frequencies represent a weighted mean across regions.
        regions is a list of region names
        weights is a list of region weights
        """

        clade_frequencies = defaultdict(dict)

        proportions = np.array(weights) / np.array(weights).sum(axis=0)

        for (region, proportion) in zip(regions, proportions):

            # Find tips within region.
            tips = [(tip.clade, tip.attr["num_date"])
                for tip in tree.get_terminals() if tip.attr["region"] == region]
            tips = np.array(sorted(tips, key=lambda row: row[1]))
            clades = tips[:, 0].astype(int)
            tip_dates = tips[:, 1].astype(float)

            # Map clade ids to their corresponding frequency matrix row index.
            clade_to_index = {clades[i]: i for i in range(len(clades))}

            # Calculate tip frequencies and normalize.
            normalized_freq_matrix_global = cls.estimate_frequencies(tip_dates, pivots, normalize_to=proportion, **kwargs)
            normalized_freq_matrix_regional = cls.estimate_frequencies(tip_dates, pivots, normalize_to=1.0, **kwargs)

            for clade in clades:
                clade_frequencies["global"][clade] = normalized_freq_matrix_global[clade_to_index[clade]]
                clade_frequencies[region][clade] = normalized_freq_matrix_regional[clade_to_index[clade]]

        for node in tree.find_clades(order="postorder"):
            if node.is_terminal():
                # Set regional frequencies for tips from different regions to zero.
                for region in regions:
                    if not node.clade in clade_frequencies[region]:
                        clade_frequencies[region][node.clade] = np.zeros_like(pivots)
            else:
                clade_frequencies["global"][node.clade] = np.array(
                    [clade_frequencies["global"][child.clade] for child in node.clades]
                ).sum(axis=0)
                for region in regions:
                    clade_frequencies[region][node.clade] = np.array(
                        [clade_frequencies[region][child.clade] for child in node.clades]
                    ).sum(axis=0)

        return clade_frequencies

if __name__=="__main__":
    plot=True
    fe = test_nested_estimator()
    #fe = test_simple_estimator()
