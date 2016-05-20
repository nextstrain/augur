import numpy as np
from H3N2 import flu
from base.prediction import tree_predictor

class flu_predictor(tree_predictor):
	"""class to train and evaluate influenza prediction models"""
	def __init__(self, *args, clade_frequencies=None, **kwarks):
		super(flu_predictor, self).__init__(*args, **kwarks)
		self.global_freqs = clade_frequencies

		self.train_intervals = [(x-2, x+2.0/12) for x in range(2006,2015)]

	def compare_dfreq_with_fixation(self):
		self.train_frequencies={}
		for tint in self.train_intervals:
			self.set_train(tint)
			self.train_frequencies[tint] = self.estimate_training_frequencies()



