# estimates clade frequencies
from __future__ import division, print_function
from scipy.interpolate import interp1d
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
            self.pivot_freq = logit_inv(x, self.pc)
            try:
                freq = interp1d(self.pivots, x, kind=self.interpolation_type,bounds_error = False, assume_sorted=True) # Throws a bug with some numpy installations, isn't necessary other than for speed.
            except:
                freq = interp1d(self.pivots, x, kind=self.interpolation_type,bounds_error = False)
            estfreq = freq(self.tps)
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
        self.pivot_freq = fix_freq(self.pivot_freq, 0.0001)

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
        self.pivot_freq[self.pivots<self.pivot_lower_cutoff] = self.fe.pivot_freq[0]
        self.pivot_freq[self.pivots>=self.pivot_upper_cutoff] = self.fe.pivot_freq[-1]


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
    def __init__(self, tree, pivots, node_filter=None, min_clades = 20, verbose=0, **kwargs):
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
        self.verbose=verbose
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

                ne = nested_frequencies(node_tps, obs_to_estimate, self.pivots, **self.kwargs)
                freq_est = ne.calc_freqs()
                for clade, tmp_freq in freq_est.iteritems():
                    if clade!="other":
                        self.frequencies[clade] = self.frequencies[node.clade]*tmp_freq
                if len(small_clades)>1:
                    for clade in small_clades:
                        if len(node.leafs):
                            frac = 1.0*len(clade.leafs)/len(node.leafs)
                        else:
                            frac = 0.0
                        self.frequencies[clade.clade] = frac*freq_est["other"]
            else:
                for clade in small_clades:
                    if len(node.leafs):
                        frac = 1.0*len(clade.leafs)/len(node.leafs)
                    else:
                        frac = 0.0
                    self.frequencies[clade.clade] = frac*self.frequencies[node.clade]


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


if __name__=="__main__":
    plot=True
    fe = test_nested_estimator()
    #fe = test_simple_estimator()
