# estimates clade frequencies
from scipy.interpolate import interp1d
import time
import numpy as np

debug = False
log_thres = 7.0

def running_average(obs, ws):
    '''
    calculates a running average
    obs     --  observations
    ws      --  window size (number of points to average)
    '''
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

def logit_transform(freq):
    return np.log(np.maximum(freq, 1e-10)/np.maximum(1e-10,(1-freq)))

def logit_inv(logit_freq):
    logit_freq[logit_freq<-log_thres]=-log_thres
    logit_freq[logit_freq>log_thres]=log_thres
    tmp_freq = np.exp(logit_freq)
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

    def __init__(self, tps, obs, pivots = None, stiffness = 20.0,
                inertia = 0.0,  tol=1e-3, pc=1e-4, ws=100, **kwargs):
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
        dt = self.tps.max()-self.tps.min()

        if pivots is None:
            self.pivot_tps = get_pivots(self.tps.min(), self.tps.max())
        elif np.isscalar(pivots):
            self.pivot_tps = np.linspace(self.tps.min()-0.01*dt, self.tps.max()+0.01*dt, pivots)
        else:
            self.pivot_tps = pivots


    def initial_guess(self, pc=0.01):
        # generate a useful initital case from a running average of the counts
        if self.ws<len(self.obs):
            tmp_vals = running_average(self.obs, self.ws)
        else:
            tmp_vals = running_average(self.obs, len(self.obs))

        tmp_interpolator = interp1d(self.tps, tmp_vals, bounds_error=False, fill_value = -1)
        pivot_freq = tmp_interpolator(self.pivot_tps)
        pivot_freq[self.pivot_tps<=tmp_interpolator.x[0]] = tmp_vals[0]
        pivot_freq[self.pivot_tps>=tmp_interpolator.x[-1]] = tmp_vals[-1]
        pivot_freq =fix_freq(pivot_freq, pc)
        return pivot_freq

    def stiffLH(self):
        freq = self.pivot_freq
        dfreq = np.diff(freq)
        dt = np.diff(self.pivot_tps)
        pq_freq = pq(fix_freq(freq,self.pc))
        # return wright fisher diffusion likelihood for frequency change.
        return -0.25*self.stiffness*(np.sum((dfreq[1:] - self.inertia*dfreq[:-1])**2/(dt[1:]*pq_freq[1:-1]))
                                    +dfreq[0]**2/(dt[0]*pq_freq[0]))


    def learn(self, initial_guess=None):
        def logLH(x):
            self.pivot_freq = x
            freq = interp1d(self.pivot_tps, self.pivot_freq, kind=self.interpolation_type)
            estfreq = fix_freq(freq(self.tps), self.pc)
            bernoulli_LH = np.sum(np.log(estfreq[self.obs])) + np.sum(np.log((1-estfreq[~self.obs])))

            stiffness_LH = self.stiffLH()
            LH = stiffness_LH + bernoulli_LH
            return -LH \
                    + 100000*(np.sum((self.pivot_freq<0)*np.abs(self.pivot_freq)) \
                    + np.sum((self.pivot_freq>1)*np.abs(self.pivot_freq-1)))

        from scipy.optimize import minimize
        if initial_guess is None:
            initial_freq = self.initial_guess(pc=0.01)
        else:
            initial_freq = initial_guess(self.pivot_tps)

        self.frequency_estimate = interp1d(self.pivot_tps, initial_freq, kind=self.interpolation_type, bounds_error=False)

        self.pivot_freq = self.frequency_estimate(self.pivot_tps)

        # determine the optimal pivot freqquencies
        self.sol = minimize(logLH, self.pivot_freq, bounds  = [(1e-5,1-1e-5) for ii in range(len(self.pivot_freq))] )
        if self.sol['success']:
            self.pivot_freq = self.sol['x']
            self.covariance = self.sol['hess_inv'].todense()
            self.conf95 = 1.96*np.sqrt(np.diag(self.covariance))
        else:
            print("Optimization failed, trying with fmin")
            #import ipdb; ipdb.set_trace()
            from scipy.optimize import fmin_powell
            self.pivot_freq = fmin_powell(logLH, self.pivot_freq, ftol = self.tol, xtol = self.tol, disp = self.verbose>0)
        self.pivot_freq = fix_freq(self.pivot_freq, 0.0001)

        # instantiate an interpolation object based on the optimal frequency pivots
        self.frequency_estimate = interp1d(self.pivot_tps, self.pivot_freq, kind=self.interpolation_type, bounds_error=False)

        if self.verbose: print "neg logLH using",len(self.pivot_tps),"pivots:", self.logLH(self.pivot_freq)

class nested_frequencies(object):
    """docstring for nested_frequencies"""
    def __init__(self, tps, obs, **kwargs):
        super(nested_frequencies, self).__init__()
        self.tps = tps
        self.obs = obs
        self.kwargs = kwargs

    def calc_freqs(self):
        sorted_obs = sorted(self.obs.items(), key=lambda x:x[1].sum(), reverse=True)

        fe = frequency_estimator(self.tps, sorted_obs[0][1], **self.kwargs)
        self.pivots = fe.pivot_tps
        self.remaining_freq = np.ones_like(self.pivots)

        self.frequencies = {}
        valid_tps = np.ones_like(self.tps, dtype=bool)
        for mut, obs in sorted_obs[:-1]:
            print(mut,'...',)
            fe = frequency_estimator(self.tps[valid_tps], obs[valid_tps], **self.kwargs)
            fe.learn()
            print('done')
            self.frequencies[mut] = self.remaining_freq * fe.pivot_freq
            self.remaining_freq *= (1.0-fe.pivot_freq)
            valid_tps = valid_tps&(~obs)

        self.frequencies[sorted_obs[-1][0]] = self.remaining_freq
        return self.frequencies

class tree_frequencies(object):
    '''
    class that estimates frequencies for nodes in the tree. each internal node is assumed
    to be named with an attribute clade, of root doesn't have such an attribute, clades
    will be numbered in preorder. Each node is assumed to have an attribute "num_date"
    '''
    def __init__(self, tree, node_filter=None, min_clades = 20, **kwargs):
        self.tree = tree
        self.min_clades = min_clades
        self.kwargs = kwargs
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
                    tps.append(node.num_date)
                    node.leafs = np.array([leaf_count])
                    leaf_count+=1
                else:
                    node.leafs = np.array([])
            else:
                node.leafs = np.concatenate([c.leafs for c in node.clades])
        self.tps = np.array(tps)

    def estimate_clade_frequencies(self):
        self.frequencies = {self.tree.root.clade:1.0}
        for node in self.tree.get_nonterminals(order='preorder'):
            print("Estimating frequencies of children of node ",node.clade)
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

                ne = nested_frequencies(node_tps, obs_to_estimate, **self.kwargs)
                freq_est = ne.calc_freqs()
                for clade, tmp_freq in freq_est.iteritems():
                    if clade!="other":
                        self.frequencies[clade] = self.frequencies[node.clade]*tmp_freq
                if len(small_clades)>1:
                    for clade in small_clades:
                        frac = 1.0*len(clade.leafs)/len(node.leafs)
                        self.frequencies[clade.clade] = frac*freq_est["other"]
            else:
                for clade in small_clades:
                    frac = 1.0*len(clade.leafs)/len(node.leafs)
                    self.frequencies[clade.clade] = frac*self.frequencies[node.clade]


class alignment_frequencies(object):

    def __init__(self, aln, tps, **kwargs):
        self.aln = np.array(aln)
        self.tps = np.array(tps)
        self.kwargs = kwargs

    def estimate_genotype_frequency(self, gt):
        match = []
        for pos, state in gt:
            match.append(aln[:,pos]==state)
        obs = match.all(axis=0)

        fe = frequency_estimator(zip(self.tps, obs), **kwargs)
        fe.learn()
        return fe.frequency_estimate


    def mutation_frequencies(self, min_freq=0.01, ignore_gap=True):
        alphabet = np.unique(self.aln)
        af = np.zeros((len(alphabet), self.aln.shape[1]))
        for ni, n in enumerate(alphabet):
            af[ni] = (self.aln==n).mean(axis=0)


        if ignore_gap:
            minor_freqs = af[alphabet!='-'].sum(axis=0) - af.max(axis=0)
        else:
            minor_freqs = 1.0 - af.max(axis=0)
        self.frequencies = {}
        for pos in np.where(minor_freqs>min_freq)[0]:
            nis = np.argsort(af[:,pos])[::-1]
            nis = nis[af[nis,pos]>0]
            muts = alphabet[nis]

            obs = {}
            for ni, mut in zip(nis, muts):
                if af[ni,pos]>min_freq and af[ni,pos]<1-min_freq:
                    obs[(pos, mut)] = self.aln[:,pos]==mut
                else:
                    break
            if len(obs)!=len(nis):
                tmp = ~np.any(obs.values(), axis=0)
                if any(tmp):
                    if len(obs)==len(nis)-1:
                        obs[(pos, muts[-1])] = tmp
                    else:
                        obs[(pos, 'other')] = tmp
            print("Estimating freuencies of position:", pos)
            print("Variants found at frequency:", [(k,o.mean()) for k,o in obs.iteritems()])
            ne = nested_frequencies(self.tps, obs, **self.kwargs)
            self.frequencies.update(ne.calc_freqs())


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



