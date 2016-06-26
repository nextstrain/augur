import numpy as np
#from flu import H3N2
from base.prediction import tree_predictor, LBI
from H3N2 import flu_process

class flu_predictor(tree_predictor):
    """class to train and evaluate influenza prediction models"""
    def __init__(self, clade_frequencies=None, pivots=None, **kwarks):
        super(flu_predictor, self).__init__(**kwarks)
        self.global_freqs = clade_frequencies
        self.global_pivots = pivots

        self.train_intervals = [(x-2, x+2.0/12) for x in range(2009,2017)]

    def compare_dfreq_with_fixation(self, **kwarks):
        self.train_frequencies={}
        for tint in self.train_intervals:
            self.set_train(tint)
            self.train_frequencies[tint] = self.estimate_training_frequencies(**kwarks)


    def calculate_LBI(self):
        self.train_frequencies={}
        for node in self.tree.find_clades():
            node.LBI={}
        for tint in self.train_intervals:
            self.set_train(tint)
            LBI(tree, tau=0.0005, attr='lbi')
            for node in self.tree.find_clades():
                node.LBI[tint] = node.lbi


    def frequency_prediction(self, train_interval):
        from collections import defaultdict
        self.predicted_frequencies = {}
        nodes_by_time = defaultdict(list)
        for node in self.tree.find_clades():
            nodes_by_time[node.numdate].append(node)
        self.all_times = sorted(nodes_by_time.keys())

        t_lower = self.all_times[0]
        self.predictions = defaultdict(dict)
        for t_cut in self.all_times[1:]:
            if t_cut<train_interval[0] or t_cut>=train_interval[1]:
                continue
            cut = []
            for node in self.tree.find_clades():
                if (node.numdate>=t_cut) and (node.up.numdate<t_cut):
                    cut.append(node)

            t0 = self.train_frequencies[train_interval][0][-2]
            f0 = {node:self.train_frequencies[train_interval][1][node.clade][-2] for node in cut}
            fitness = {node:node.fitness[train_interval] for node in cut}

            tmp_freqs = {node:np.exp((self.global_pivots-t0)*node.fitness[train_interval])*f0[node] for node in cut}
            total_freq = np.sum(tmp_freqs.values(), axis=0)
            for node in cut:
                self.predictions[t_cut][node] = tmp_freqs[node]/total_freq

    def plot_prediction(self,t_cut):
        ti = np.searchsorted(self.all_times, t_cut)
        t_cut_val = self.all_times[ti]
        for node in self.predictions:
            if t_cut_val in self.predictions[node]:
                plt.plot(self.global_pivots, self.predictions[node][t_cut_val])
                plt.plot(self.global_pivots, self.global_freqs[node.clade])


    def slope_fitness(self):
        for node in self.tree.find_clades():
            node.fitness={}
        for train_interval,(pivots, freq) in self.train_frequencies.iteritems():
            dt = pivots[-1]-pivots[-2]
            for node in self.tree.find_clades():
                slope = (freq[node.clade][-1] - freq[node.clade][-2])/dt
                node.fitness[train_interval] = slope



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

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('--load', action='store_true', help = 'recover from file')

    params = parser.parse_args()
    fname = sorted(glob.glob('../nextstrain-db/data/flu_h3n2*fasta'))[-1]
    titer_fname = sorted(glob.glob('../nextstrain-db/data/h3n2*text'))[-1]

    flu = flu_process(method='SLSQP', dtps=2.0, stiffness=20, inertia=0.9,
                      fname = fname)
    if params.load:
        flu.load()
        flu.estimate_tree_frequencies()
    else:
        pass

    flu.HI_model(titer_fname)


    flu_pred = flu_predictor(tree=flu.tree.tree, seqs=flu.seqs,
                             clade_frequencies = flu.tree_frequencies,
                             pivots=flu.tree_pivots,inertia=1.0, stiffness=20)

    flu_pred.compare_dfreq_with_fixation()

    flu_pred.plot_all()
    global_freqs, interval_freqs = flu_pred.delta_freq_vs_success()

    failed = global_freqs[:,0]<0.1
    fixed = global_freqs[:,0]>0.9
    failed_peaks = global_freqs[failed,1]
    plt.figure()
    plt.plot(sorted(failed_peaks), np.linspace(1,0,len(failed_peaks)), label='failed')
    plt.plot(sorted(global_freqs[:,1]), np.linspace(1,0,len(global_freqs)), label='all')
    plt.xlabel('peak frequency')
    plt.ylabel('fraction with higher peak')

    failed_jumps = global_freqs[failed,2]
    fixed_jumps = global_freqs[fixed,2]
    plt.figure()
    plt.plot(sorted(failed_jumps), np.linspace(1,0,len(failed_jumps)), label='failed')
    plt.plot(sorted(fixed_jumps), np.linspace(1,0,len(fixed_jumps)), label='fixed')
    plt.xlabel('max jump')
    plt.ylabel('fraction with larger jump')






