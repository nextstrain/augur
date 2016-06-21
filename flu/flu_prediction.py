import numpy as np
#from flu import H3N2
from base.prediction import tree_predictor
from base.titer_model import tree_model
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


if __name__=="__main__":
    import argparse, glob

    parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
    parser.add_argument('--load', action='store_true', help = 'recover from file')

    params = parser.parse_args()
    fname = sorted(glob.glob('../nextstrain-db/data/flu_h3n2*fasta'))[-1]

    flu = flu_process(method='SLSQP', dtps=2.0, stiffness=20, inertia=0.9,
                      fname = fname)
    if params.load:
        flu.load()
    else:


    flu.estimate_tree_frequencies()
    titer_fname = sorted(glob.glob('../nextstrain-db/data/h3n2*text'))[-1]
    flu_titers = tree_model(flu.tree.tree, titer_fname = titer_fname)
    flu_titers.prepare()
    flu_titers.train()

    flu_pred = flu_predictor(tree=flu.tree.tree, seqs=flu.seqs,
                             clade_frequencies = flu.tree_frequencies,
                             pivots=flu.tree_pivots,inertia=1.0, stiffness=20)

    flu_pred.compare_dfreq_with_fixation()

    flu_pred.plot_all()
