import numpy as np
#from flu import H3N2
from base.prediction import tree_predictor
from base.titer_model import tree_model

class flu_predictor(tree_predictor):
    """class to train and evaluate influenza prediction models"""
    def __init__(self, clade_frequencies=None, **kwarks):
        super(flu_predictor, self).__init__(**kwarks)
        self.global_freqs = clade_frequencies

        self.train_intervals = [(x-2, x+2.0/12) for x in range(2009,2017)]

    def compare_dfreq_with_fixation(self):
        self.train_frequencies={}
        for tint in self.train_intervals:
            self.set_train(tint)
            self.train_frequencies[tint] = self.estimate_training_frequencies()




if __name__=="__main__":
    #   flu.estimate_tree_frequencies()
    flu_titers = tree_model(flu.tree.tree, titer_fname = '../../nextflu2/data/H3N2_HI_titers.txt')
    flu_titers.prepare()
    flu_titers.train()

    flu_pred = flu_predictor(tree=flu.tree.tree, seqs=flu.seqs,
                             clade_frequencies = flu.tree_frequencies)
