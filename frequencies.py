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

    def __init__(self, observations, pivots = None, stiffness = 20.0,
                inertia = 0.0,  tol=1e-3, pc=1e-4, **kwargs):
        self.tps = np.array([x[0] for x in observations])
        self.obs = np.array([x[1]>0 for x in observations])
        self.stiffness = stiffness
        self.inertia = inertia
        self.interpolation_type = 'linear'
        self.tol = tol
        self.reg = 1e-6
        self.pc = pc
        self.verbose = 0

        if pivots is None:
            self.pivot_tps = get_pivots(self.tps[0], self.tps[-1])
        elif np.isscalar(pivots):
            self.pivot_tps = np.linspace(self.tps[0], self.tps[-1], pivots)
        else:
            self.pivot_tps = pivots

    def initial_guess(self, ws=50):
        # generate a useful initital case from a running average of the counts
        tmp_vals = running_average(self.obs, ws)
        tmp_interpolator = interp1d(self.tps, tmp_vals, bounds_error=False, fill_value = -1)
        pivot_freq = tmp_interpolator(self.pivot_tps)
        pivot_freq[self.pivot_tps<=tmp_interpolator.x[0]] = tmp_vals[0]
        pivot_freq[self.pivot_tps>=tmp_interpolator.x[-1]] = tmp_vals[-1]
        pivot_freq =fix_freq(pivot_freq, 0.95)
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
            return -LH/len(self.obs) \
                    + 100000*(np.sum((self.pivot_freq<0)*np.abs(self.pivot_freq)) \
                    + np.sum((self.pivot_freq>1)*np.abs(self.pivot_freq-1)))

        from scipy.optimize import minimize
        if initial_guess is None:
            initial_freq = fix_freq(self.initial_guess(ws=2*(min(50,len(self.obs))//2)), 0.01)
        else:
            initial_freq = initial_guess(self.pivot_tps)

        self.frequency_estimate = interp1d(self.pivot_tps, initial_freq, kind=self.interpolation_type, bounds_error=False)

        self.pivot_freq = self.frequency_estimate(self.pivot_tps)

        # determine the optimal pivot freqquencies
        self.sol = minimize(logLH, self.pivot_freq, bounds  = [(0,1) for ii in range(len(self.pivot_freq))] )
        if self.sol['success']:
            self.pivot_freq = self.sol['x']
        else:
            print("Optimization failed, trying with fmin")
            from scipy.optimize import fmin_powell
            self.pivot_freq = fmin_powell(logLH, self.pivot_freq, ftol = self.tol, xtol = self.tol, disp = self.verbose>0)
        self.pivot_freq = fix_freq(self.pivot_freq, 0.0001)

        # instantiate an interpolation object based on the optimal frequency pivots
        self.frequency_estimate = interp1d(self.pivot_tps, self.pivot_freq, kind=self.interpolation_type, bounds_error=False)

        if self.verbose: print "neg logLH using",len(self.pivot_tps),"pivots:", self.logLH(self.pivot_freq)


class tree_frequencies(object):
    def __init__(self, tree):
        self.tree = tree

class alignment_frequencies(object):

    def __init__(self, aln):
        self.aln = aln

def test():
    import matplotlib.pyplot as plt
    tps = np.sort(100 * np.random.uniform(size=100))
    freq_traj = [0.1]
    stiffness=100
    s=-0.02
    for dt in np.diff(tps):
        freq_traj.append(freq_traj[-1]*np.exp(-s*dt)+np.sqrt(2*np.max(0,freq_traj[-1]*(1-freq_traj[-1]))*dt/stiffness)*np.random.normal())
    obs = np.random.uniform(size=tps.shape)<freq_traj

    fe = frequency_estimator(zip(tps, obs), pivots=20, stiffness=stiffness)
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

if __name__=="__main__":
    plot=True
    fe = test()



