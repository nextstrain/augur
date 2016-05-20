from __future__ import division, print_function
import numpy as np
import time
from collections import defaultdict
from base.io_util import myopen
from itertools import izip
import pandas as pd

def LBI(tree, tau=0.001, attr='lbi'):
'''
traverses the tree in postorder and preorder to calculate the
up and downstream tree length exponentially weighted by distance.
then adds them as LBI
tree -- Biopython tree
tau  -- time scale for tree length integration
attr -- the attribute name used to store the result
'''
# traverse the tree in postorder (children first) to calculate msg to parents
for node in tree.find_clades(order='postorder'):
    node.down_polarizer = 0
    node.up_polarizer = 0
    for child in node.clades():
        node.up_polarizer += child.up_polarizer
    bl =  node.branch_length/tau
    node.up_polarizer *= np.exp(-bl)
    if node.train: node.up_polarizer += tau*(1-np.exp(-bl))

# traverse the tree in preorder (parents first) to calculate msg to children
for node in tree.get_nonterminals(order='preorder'):
    for child1 in node.child_nodes():
        child1.down_polarizer = node.down_polarizer
        for child2 in node.child_nodes():
            if child1!=child2:
                child1.down_polarizer += child2.up_polarizer

        bl =  child1.edge_length/tau
        child1.down_polarizer *= np.exp(-bl)
        if child1.train: child1.down_polarizer += tau*(1-np.exp(-bl))

# go over all nodes and calculate the LBI (can be done in any order)
for node in tree.find_clades():
    tmp_LBI = node.down_polarizer
    for child in node.child_nodes():
        tmp_LBI += child.up_polarizer
    node.__setattr__(attr, transform(tmp_LBI))



class tree_predictor(object):
    """class implementing basic methods to extrapolate genetic diversity
    into the future. specific predictive features need to be implemented
    by the subclasses"""
    def __init__(self, tree, seqs, features=[], **kwarks):
        super(tree_predictor, self).__init__()
        self.tree = tree
        self.seqs = seqs
        self.features = features


    def set_train(train_interval, train_filter=None):
        '''
        mark all nodes in tree as test, train, or neither
        train_interval -- (start_date,  stop_date) as numerical date
        '''
        if train_filter is None: train_filter = lambda x:True
        self.train_interval = train_interval
        in_interval = lambda x,y: (x>=y[0])&(x<y[1])
        for node in self.tree.find_clades():
            node.train = in_interval(node.num_date, train_interval[0])&train_filter(node)


    def set_test(test_interval, test_filter=None):
        '''
        mark all nodes in tree as test, train, or neither
        test_interval  -- (start_date,  stop_date) as numerical date
        '''
        if test_filter is None: test_filter = lambda x:True
        self.test_interval  = test_interval
        in_interval = lambda x,y: (x>=y[0])&(x<y[1])
        for node in self.tree.find_clades():
            node.test = in_interval(node.num_date, test_interval[0])&test_filter(node)


    def estimate_training_frequencies(self):
        from frequencies import tree_frequencies
        npivots = int((self.train_interval[1]-self.test_interval[0])*12)
        pivots=np.linspace(*self.train_interval, 12)
        fe = tree_frequencies(self.tree, pivots, node_filter=lambda x:x.train, min_clades=10)
        # dictionary containing frequencies of all clades.
        # the keys are the node.clade attributes
        return fe.frequencies

