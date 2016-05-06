# clean, reroot, ladderize newick tree
# output to tree.json
from __future__ import division, print_function
import numpy as np
import time, os, gzip
from collections import defaultdict
from nextstrain.io_util import myopen
from matplotlib import pyplot as plt
from itertools import izip
from H3N2 import fix_name
import pandas as pd


class HI_tree(object):

    def __init__(self, tree,  HI_fname = 'data/HI_titers.txt',serum_Kc = 0.0,
                 min_aa_muts = 0,**kwargs):
        self.HI_fname = HI_fname
        self.tree = tree
        self.node_lookup = {n.name:n for n in tree.get_terminals()}
        self.node_lookup.update({n.name.upper():n for n in tree.get_terminals()})
        self.node_lookup.update({n.name.lower():n for n in tree.get_terminals()})
        self.tree.root.parent_node=None
        for node in self.tree.get_nonterminals():
            for c in node.clades:
                c.parent_node = node

        if "excluded_tables" in kwargs:
            self.excluded_tables = kwargs["excluded_tables"]
        else:
            self.excluded_tables = []
        self.HI, tmp, sources = self.read_HI_titers(HI_fname)
        self.sources = list(sources)
        self.serum_Kc = serum_Kc
        self.tree_graph = None
        self.min_aa_muts = min_aa_muts
        self.serum_potency = {}
        self.virus_effect = {}

    def read_HI_titers(self, fname):
        strains = set()
        measurements = defaultdict(list)
        sources = set()
        with myopen(fname, 'r') as infile:
            for line in infile:
                entries = line.strip().split()
                test, ref_virus, serum, src_id, val = (entries[0], entries[1],entries[2],
                                                        entries[3], float(entries[4]))
                ref = (ref_virus, serum)
                if src_id not in self.excluded_tables:
                    try:
                        measurements[(test, (ref_virus, serum))].append(val)
                        strains.update([test, ref])
                        sources.add(src_id)
                    except:
                        print(line.strip())
        return measurements, strains, sources

    def normalize(self, ref, val):
        consensus_func = np.mean
        return consensus_func(np.log2(self.autologous_titers[ref]['val'])) \
                - consensus_func(np.log2(val))

    def normalize_HI(self):
        '''
        convert the HI measurements into the log2 difference between the average
        HI titer measured between test virus and reference serum and the average
        homologous titer. all measurements relative to sera without homologous titer
        are excluded
        '''
        self.HI_normalized = {}
        self.HI_raw = {}
        self.measurements_per_serum = defaultdict(int)
        sera = set()
        ref_strains = set()
        HI_strains = set()
        all_per_serum = defaultdict(list)
        autologous = defaultdict(list)
        for (test, ref), val in self.HI.iteritems():
            if test.upper() in self.node_lookup and ref[0].upper() in self.node_lookup:
                HI_strains.add(test.upper())
                HI_strains.add(ref[0].upper())
                all_per_serum[ref].append(val)
                if ref[0]==test:
                    autologous[ref].append(val)

        self.autologous_titers = {}
        for serum in all_per_serum:
            if serum in autologous:
                self.autologous_titers[serum] = {'val':autologous[serum], 'autologous':True}
                print("autologous titer found for",serum)
            else:
                if len(all_per_serum[serum])>10:
                    self.autologous_titers[serum] = {'val':np.max(all_per_serum[serum]), 'autologous':False}
                    print(serum,": using max titer instead of autologous,", np.max(all_per_serum[serum]))
                else:
                    print("discarding",serum,"since there are only ",len(all_per_serum[serum]),'measurements')

        for (test, ref), val in self.HI.iteritems():
            if test.upper() in self.node_lookup and ref[0].upper() in self.node_lookup:
                if ref in self.autologous_titers:
                    sera.add(ref)
                    ref_strains.add(ref[0])
                    self.HI_normalized[(test, ref)] = self.normalize(ref, val)
                    self.HI_raw[(test, ref)] = np.median(val)
                    self.measurements_per_serum[ref]+=1
                else:
                    pass
                    #print "no homologous titer found:", ref

        self.sera = list(sera)
        self.ref_strains = list(ref_strains)
        self.HI_strains = list(HI_strains)

    def add_mutations(self):
        '''
        add amino acid mutation count to the tree
        '''
        self.tree.root.n_aa_muts= 0
        for node in self.tree.find_clades(order='postorder'):
            if node is not self.tree.root:
                node.n_aa_muts = 0
                for prot in node.aa_muts:
                    if node.aa_muts[prot]!='':
                        node.n_aa_muts += len(node.aa_muts[prot].split(','))

    def mark_HI_splits(self):
        # flag all branches on the tree with HI_strain = True if they lead to strain with titer data
        for leaf in self.tree.get_terminals():
            if leaf.strain.upper() in self.HI_strains:
                leaf.serum = leaf.strain.upper() in self.ref_strains
                leaf.HI_info = 1
            else:
                leaf.serum, leaf.HI_info=False, 0

        for node in self.tree.get_nonterminals(order='postorder'):
            node.HI_info = sum([c.HI_info for c in node.clades])
            node.serum= False

        # combine sets of branches that span identical sets of HI measurements
        self.HI_split_count = 0  # HI measurment split counter
        self.HI_split_to_branch = defaultdict(list)
        for node in self.tree.find_clades(order='preorder'):
            if self.map_to_tree:
                node.dHI, node.cHI, node.mHI, node.constraints = 0, 0, 0, 0
            if node.HI_info>1:
                node.HI_branch_index = self.HI_split_count
                self.HI_split_to_branch[node.HI_branch_index].append(node)
                # at a bi- or multifurcation, increase the split count and HI index
                # either individual child branches have enough HI info be counted,
                # or the pre-order node iteraction will move towards the root
                if sum([c.HI_info>0 for c in node.clades])>1:
                    self.HI_split_count+=1
                elif node.is_terminal():
                    self.HI_split_count+=1

        self.genetic_params = self.HI_split_count
        print ("# of reference strains:",len(self.sera), "# of branches with HI constraint", self.HI_split_count)

    def get_path_no_terminals(self, v1, v2):
        '''
        returns the path between two tips in the tree excluding the terminal branches.
        '''
        if v1 in self.node_lookup and v2 in self.node_lookup:
            p1 = [self.node_lookup[v1]]
            p2 = [self.node_lookup[v2]]
            for tmp_p in [p1,p2]:
                while tmp_p[-1].parent_node != self.tree.root:
                    tmp_p.append(tmp_p[-1].parent_node)
                tmp_p.append(self.tree.root)
                tmp_p.reverse()

            for pi, (tmp_v1, tmp_v2) in enumerate(izip(p1,p2)):
                if tmp_v1!=tmp_v2:
                    break
            path = [n for n in p1[pi:] if n.HI_info>1] + [n for n in p2[pi:] if n.HI_info>1]
        else:
            path = None
        return path

    def get_mutations(self, strain1, strain2):
        ''' return amino acid mutations between viruses specified by strain names as tuples (HA1, F159S) '''
        if strain1 in self.node_lookup and strain2 in self.node_lookup:
            node1 = self.node_lookup[strain1].parent_node
            node2 = self.node_lookup[strain2].parent_node
            if node1 is None or node2 is None:
                return None
            else:
                return self.get_mutations_nodes(node1, node2)
        else:
            return None

    def get_mutations_nodes(self, node1, node2):
        muts = []
        for prot in node1.aa_seq:
            seq1 = node1.aa_seq[prot]
            seq2 = node2.aa_seq[prot]
            muts.extend([(prot, aa1+str(pos+1)+aa2) for pos, (aa1, aa2) in enumerate(izip(seq1, seq2)) if aa1!=aa2])
        return muts

    def make_seqgraph(self, colin_thres = 5):
        '''
        code amino acid differences between sequences into a matrix
        the matrix has dimensions #measurements x #observed mutations
        '''
        seq_graph = []
        HI_dist = []
        weights = []
        # count how often each mutation separates a reference test virus pair
        self.mutation_counter = defaultdict(int)
        for (test, ref), val in self.train_HI.iteritems():
            muts = self.get_mutations(ref[0], test)
            if muts is None:
                continue
            for mut in muts:
                self.mutation_counter[mut]+=1

        # make a list of mutations deemed relevant via frequency thresholds
        relevant_muts = []
        min_count= 10
        min_freq = 1.0*min_count/len(self.viruses)
        for mut, count in self.mutation_counter.iteritems():
            gene = mut[0]
            pos = int(mut[1][1:-1])-1
            aa1, aa2 = mut[1][0],mut[1][-1]
            if gene=='HA1' and count>min_count and \
                self.aa_frequencies[gene][self.aa_alphabet.index(aa1),pos]>min_freq and\
                self.aa_frequencies[gene][self.aa_alphabet.index(aa2),pos]>min_freq:
                relevant_muts.append(mut)

        relevant_muts.sort(key = lambda x:int(x[1][1:-1]))
        self.relevant_muts = relevant_muts
        self.genetic_params = len(relevant_muts)
        n_params = self.genetic_params + len(self.sera) + len(self.HI_strains)

        # loop over all measurements and encode the HI model as [0,1,0,1,0,0..] vector:
        # 1-> mutation present, 0 not present, same for serum and virus effects
        for (test, ref), val in self.train_HI.iteritems():
            if not np.isnan(val):
                try:
                    muts = self.get_mutations(ref[0], test)
                    if muts is None:
                        continue
                    tmp = np.zeros(n_params) # zero vector, ones will be filled in
                    # determine branch indices on path
                    mutation_indices = np.unique([relevant_muts.index(mut) for mut in muts
                                         if mut in relevant_muts])
                    if len(mutation_indices): tmp[mutation_indices] = 1
                    # add serum effect for heterologous viruses
                    if test!=ref[0]:
                        tmp[len(relevant_muts)+self.sera.index(ref)] = 1
                    # add virus effect
                    tmp[len(relevant_muts)+len(self.sera)+self.HI_strains.index(test)] = 1
                    # append model and fit value to lists seq_graph and HI_dist
                    seq_graph.append(tmp)
                    HI_dist.append(val)
                    # for each measurment (row in the big matrix), attach weight that accounts for representation of serum
                    weights.append(1.0/(1.0 + self.serum_Kc*self.measurements_per_serum[ref]))
                except:
                    import pdb; pdb.set_trace()
                    print(test, ref, "ERROR")

        # convert to numpy arrays and save product of tree graph with its transpose for future use
        self.weights = np.sqrt(weights)
        self.HI_dist =  np.array(HI_dist)*self.weights
        self.tree_graph = (np.array(seq_graph).T*self.weights).T
        if colin_thres is not None:
            self.collapse_colinear_mutations(colin_thres)
        self.TgT = np.dot(self.tree_graph.T, self.tree_graph)
        print ("Found", self.tree_graph.shape, "measurements x parameters")

    def collapse_colinear_mutations(self, colin_thres):
        '''
        find colinear columns of the design matrix, collapse them into clusters
        '''
        TT = self.tree_graph[:,:self.genetic_params].T
        mutation_clusters = []
        n_measurements = self.tree_graph.shape[0]
        # a greedy algorithm: if column is similar to existing cluster -> merge with cluster, else -> new cluster
        for col, mut in izip(TT, self.relevant_muts):
            col_found = False
            for cluster in mutation_clusters:
                # similarity is defined as number of measurements at whcih the cluster and column differ
                if np.sum(col==cluster[0])>=n_measurements-colin_thres:
                    cluster[1].append(mut)
                    col_found=True
                    print("adding",mut,"to cluster ",cluster[1])
                    break
            if not col_found:
                mutation_clusters.append([col, [mut]])
        print("dimensions of old design matrix",self.tree_graph.shape)
        self.tree_graph = np.hstack((np.array([c[0] for c in mutation_clusters]).T,
                                     self.tree_graph[:,self.genetic_params:]))
        self.genetic_params = len(mutation_clusters)
        # use the first mutation of a cluster to index the effect
        # make a dictionary that maps this effect to the cluster
        self.mutation_clusters = {c[1][0]:c[1] for c in mutation_clusters}
        self.relevant_muts = [c[1][0] for c in mutation_clusters]
        print("dimensions of new design matrix",self.tree_graph.shape)

    def make_treegraph(self):
        '''
        code the path between serum and test virus of each HI measurement into a matrix
        the matrix has dimensions #measurements x #tree branches with HI info
        if the path between test and serum goes through a branch, the corresponding matrix element is 1, 0 otherwise
        '''
        tree_graph = []
        HI_dist = []
        weights = []
        # mark HI splits have to have been run before, assigning self.HI_split_count
        n_params = self.HI_split_count + len(self.sera) + len(self.HI_strains)
        self.genetic_params = self.HI_split_count
        for (test, ref), val in self.train_HI.iteritems():
            if not np.isnan(val):
                if True: #try:
                    if ref[0] in self.node_lookup and test in self.node_lookup\
                        and self.node_lookup[ref[0]].parent_node is not None\
                        and self.node_lookup[test].parent_node is not None:
                        path = self.get_path_no_terminals(test, ref[0])
                        tmp = np.zeros(n_params)
                        # determine branch indices on path
                        if type(self.min_aa_muts)==int:
                            branches = np.unique([c.HI_branch_index for c in path
                                                 if (hasattr(c, 'HI_branch_index') and
                                                     c.n_aa_muts>=self.min_aa_muts)])
                        elif self.min_aa_muts=='epi':
                            branches = np.unique([c.HI_branch_index for c in path
                                                 if (hasattr(c, 'HI_branch_index') and c.parent_node.ep<c.ep)])
                        elif self.min_aa_muts=='rbs':
                            branches = np.unique([c.HI_branch_index for c in path
                                                 if (hasattr(c, 'HI_branch_index') and c.parent_node.rb<c.rb)])
                        else:
                            branches = np.unique([c.HI_branch_index for c in path
                                                 if hasattr(c, 'HI_branch_index') ])
                        if len(branches): tmp[branches] = 1
                        # add serum effect for heterologous viruses
                        if ref[0]!=test:
                            tmp[self.HI_split_count+self.sera.index(ref)] = 1
                        # add virus effect
                        tmp[self.HI_split_count+len(self.sera)+self.HI_strains.index(test)] = 1
                        # append model and fit value to lists tree_graph and HI_dist
                        tree_graph.append(tmp)
                        HI_dist.append(val)
                        weights.append(1.0/(1.0 + self.serum_Kc*self.measurements_per_serum[ref]))
                #except:
                #    import ipdb; ipdb.set_trace()
                #    print test, ref, "ERROR"

        # convert to numpy arrays and save product of tree graph with its transpose for future use
        self.weights = np.sqrt(weights)
        self.HI_dist =  np.array(HI_dist)*self.weights
        self.tree_graph = (np.array(tree_graph).T*self.weights).T
        self.TgT = np.dot(self.tree_graph.T, self.tree_graph)
        print ("Found", self.tree_graph.shape, "measurements x parameters")

    def fit_func(self):
        return np.mean( (self.HI_dist - np.dot(self.tree_graph, self.params))**2 )

    def fit_l1reg(self):
        from cvxopt import matrix, solvers
        n_params = self.tree_graph.shape[1]
        HI_sc = self.genetic_params if self.map_to_tree else len(self.relevant_muts)
        n_sera = len(self.sera)
        n_v = len(self.HI_strains)

        # set up the quadratic matrix containing the deviation term (linear xterm below)
        # and the l2-regulatization of the avidities and potencies
        P1 = np.zeros((n_params+HI_sc,n_params+HI_sc))
        P1[:n_params, :n_params] = self.TgT
        for ii in xrange(HI_sc, HI_sc+n_sera):
            P1[ii,ii]+=self.lam_pot
        for ii in xrange(HI_sc+n_sera, n_params):
            P1[ii,ii]+=self.lam_avi
        P = matrix(P1)

        # set up cost for auxillary parameter and the linear cross-term
        q1 = np.zeros(n_params+HI_sc)
        q1[:n_params] = -np.dot( self.HI_dist, self.tree_graph)
        q1[n_params:] = self.lam_HI
        q = matrix(q1)

        # set up linear constraint matrix to regularize the HI parametesr
        h = matrix(np.zeros(2*HI_sc))   # Gw <=h
        G1 = np.zeros((2*HI_sc,n_params+HI_sc))
        G1[:HI_sc, :HI_sc] = -np.eye(HI_sc)
        G1[:HI_sc:, n_params:] = -np.eye(HI_sc)
        G1[HI_sc:, :HI_sc] = np.eye(HI_sc)
        G1[HI_sc:, n_params:] = -np.eye(HI_sc)
        G = matrix(G1)
        W = solvers.qp(P,q,G,h)
        self.params = np.array([x for x in W['x']])[:n_params]
        print ("rms deviation prior to relax=",np.sqrt(self.fit_func()))
        return self.params

    def fit_nnls(self):
        from scipy.optimize import nnls
        return nnls(self.tree_graph, self.HI_dist)[0]

    def fit_nnl2reg(self):
        from cvxopt import matrix, solvers
        n_params = self.tree_graph.shape[1]
        P = matrix(np.dot(self.tree_graph.T, self.tree_graph) + self.lam_HI*np.eye(n_params))
        q = matrix( -np.dot( self.HI_dist, self.tree_graph))
        h = matrix(np.zeros(n_params)) # Gw <=h
        G = matrix(-np.eye(n_params))
        W = solvers.qp(P,q,G,h)
        return np.array([x for x in W['x']])

    def fit_nnl1reg(self):
        from cvxopt import matrix, solvers
        n_params = self.tree_graph.shape[1]
        HI_sc = self.genetic_params
        n_sera = len(self.sera)
        n_v = len(self.HI_strains)

        # set up the quadratic matrix containing the deviation term (linear xterm below)
        # and the l2-regulatization of the avidities and potencies
        P1 = np.zeros((n_params,n_params))
        P1[:n_params, :n_params] = self.TgT
        for ii in xrange(HI_sc, HI_sc+n_sera):
            P1[ii,ii]+=self.lam_pot
        for ii in xrange(HI_sc+n_sera, n_params):
            P1[ii,ii]+=self.lam_avi
        P = matrix(P1)

        # set up cost for auxillary parameter and the linear cross-term
        q1 = np.zeros(n_params)
        q1[:n_params] = -np.dot(self.HI_dist, self.tree_graph)
        q1[:HI_sc] += self.lam_HI
        q = matrix(q1)

        # set up linear constraint matrix to enforce positivity of the
        # dHIs and bounding of dHI by the auxillary parameter
        h = matrix(np.zeros(HI_sc))     # Gw <=h
        G1 = np.zeros((HI_sc,n_params))
        G1[:HI_sc, :HI_sc] = -np.eye(HI_sc)
        G = matrix(G1)
        W = solvers.qp(P,q,G,h)
        self.params = np.array([x for x in W['x']])[:n_params]
        print ("rms deviation prior to relax=",np.sqrt(self.fit_func()))
        # redo the linear cost relaxing terms that seem to be relevant to avoid
        # compression of the fit. 0.2 seems to be a good cut-off, linear tune to zero
        #q1[n_params:] = self.lam_HI*(1-5.0*np.minimum(0.2,sol[:HI_sc]))
        #q = matrix(q1)
        #W = solvers.qp(P,q,G,h)
        #sol = np.array([x for x in W['x']])[:n_params]
        #self.params=sol
        #print "rms deviation after relax=",np.sqrt(self.fit_func())
        return self.params

    def prepare_HI_map(self):
        '''
        normalize the HI measurements, split the data into training and test sets
        and determine which branches on the tree are transversed by HI measurements
        '''
        from random import sample
        self.normalize_HI()
        self.add_mutations()
        self.mark_HI_splits()
        if self.training_fraction<1.0: # validation mode, set aside a fraction of measurements to validate the fit
            self.test_HI, self.train_HI = {}, {}
            if self.subset_strains:    # exclude a fraction of test viruses
                tmp = set(self.HI_strains)
                tmp.difference_update(self.ref_strains) # don't use references viruses in the set to sample from
                training_strains = sample(tmp, int(self.training_fraction*len(tmp)))
                for tmpstrain in self.ref_strains:      # add all reference viruses to the training set
                    if tmpstrain not in training_strains:
                        training_strains.append(tmpstrain)
                for key, val in self.HI_normalized.iteritems():
                    if key[0] in training_strains:
                        self.train_HI[key]=val
                    else:
                        self.test_HI[key]=val
            else: # simply use a fraction of all measurements for testing
                for key, val in self.HI_normalized.iteritems():
                    if np.random.uniform()>self.training_fraction:
                        self.test_HI[key]=val
                    else:
                        self.train_HI[key]=val
        else: # without the need for a test data set, use the entire data set for training
            self.train_HI = self.HI_normalized

        # if data is to censored by date, subset the data set and reassign sera, reference strains, and test viruses
        if self.cutoff_date is not None:
            prev_years = 6 # number of years prior to cut-off to use when fitting date censored data
            self.train_HI = {key:val for key,val in self.train_HI.iteritems()
                            if self.node_lookup[key[0]].num_date<=self.cutoff_date and
                               self.node_lookup[key[1][0]].num_date<=self.cutoff_date and
                               self.node_lookup[key[0]].num_date>self.cutoff_date-prev_years and
                               self.node_lookup[key[1][0]].num_date>self.cutoff_date-prev_years}
            sera = set()
            ref_strains = set()
            HI_strains = set()

            for test,ref in self.train_HI:
                if test.upper() in self.node_lookup and ref[0].upper() in self.node_lookup:
                    HI_strains.add(test)
                    HI_strains.add(ref[0])
                    sera.add(ref)
                    ref_strains.add(ref[0])

            self.sera = list(sera)
            self.ref_strains = list(ref_strains)
            self.HI_strains = list(HI_strains)

        # construct the design matrix depending on the model
        if self.map_to_tree:
            self.make_treegraph()
        else:
            self.make_seqgraph()

    def map_HI(self, training_fraction = 1.0, method = 'nnls', lam_HI=1.0, map_to_tree = True,
            lam_pot = 0.5, lam_avi = 3.0, cutoff_date = None, subset_strains = False, force_redo = False):
        self.map_to_tree = map_to_tree
        self.training_fraction = training_fraction
        self.subset_strains=subset_strains
        self.lam_pot = lam_pot
        self.lam_avi = lam_avi
        self.lam_HI = lam_HI
        self.cutoff_date = cutoff_date
        if self.tree_graph is None or force_redo:
            self.prepare_HI_map()

        if method=='l1reg':  # l1 regularized fit, no constraint on sign of effect
            self.params = self.fit_l1reg()
        elif method=='nnls':  # non-negative least square, not regularized
            self.params = self.fit_nnls()
        elif method=='nnl2reg': # non-negative L2 norm regularized fit
            self.params = self.fit_nnl2reg()
        elif method=='nnl1reg':  # non-negative fit, branch terms L1 regularized, avidity terms L2 regularized
            self.params = self.fit_nnl1reg()

        self.fit_error = np.sqrt(self.fit_func())
        print("method",method, "regularized by", self.lam_HI, "rms deviation=", self.fit_error)
        # for each set of branches with HI constraints, pick the branch with most aa mutations
        # and assign the dHI to that one, record the number of constraints
        if self.map_to_tree:
            for node in self.tree.find_clades(order='postorder'):
                node.dHI=0
            for HI_split, branches in self.HI_split_to_branch.iteritems():
                likely_branch = branches[np.argmax([b.n_aa_muts for b in branches])]
                likely_branch.dHI = self.params[HI_split]
                likely_branch.constraints = self.tree_graph[:,HI_split].sum()

            # integrate the tree model HI change dHI into a cumulative antigentic evolution score cHI
            for node in self.tree.find_clades(order='preorder'):
                if node.parent_node is not None:
                    node.cHI = node.parent_node.cHI + node.dHI
                else:
                    node.cHI=0
        else:
            self.mutation_effects={}
            for mi, mut in enumerate(self.relevant_muts):
                self.mutation_effects[mut] = self.params[mi]
            # integrate the mutation model change into a cumulative antigentic evolution score mHI
            for node in self.tree.find_clades(order='preorder'):
                if node.parent_node is not None:
                    node.mHI = node.parent_node.mHI \
                    + sum([0]+[self.mutation_effects[('HA1',str(mut))] for mut in node.aa_muts['HA1'].split(',')
                           if ('HA1',str(mut)) in self.mutation_effects])
                else:
                    node.mHI=0

        self.serum_potency['tree' if self.map_to_tree else 'mutation'] =\
                    {serum:self.params[self.genetic_params+ii]
                      for ii, serum in enumerate(self.sera)}
        self.virus_effect['tree' if self.map_to_tree else 'mutation'] = \
                    {strain:self.params[self.genetic_params+len(self.sera)+ii]
                  for ii, strain in enumerate(self.HI_strains)}

    def generate_validation_figures(self, method = 'nnl1reg'):
        import matplotlib.pyplot as plt

        lam_pot = self.lam_pot
        lam_avi = self.lam_avi
        lam_HI =  self.lam_HI
        # summary figures using previously determined models
        for map_to_tree, model in [(True, 'tree'), (False,'mutation')]:
            try:
                # figures showing the histograms of virus and serum effects
                plt.figure(figsize=(1.3*figheight,figheight))
                ax = plt.subplot(121)
                plt.hist(self.virus_effect[model].values(), bins=np.linspace(-2,2,21), normed=True)
                plt.xlabel('avidity', fontsize=fs)
                plt.text(0.05, 0.93,  ('tree model' if model=='tree' else 'mutation model'),
                         weight='bold', fontsize=fs, transform=plt.gca().transAxes)
                ax.set_xticks([-2,-1,0,1,2])
                ax.tick_params(axis='both', labelsize=fs)
                ax = plt.subplot(122)
                plt.hist(self.serum_potency[model].values(), bins=10, normed=True)
                plt.xlabel('potency', fontsize=fs)
                plt.tight_layout()
                ax.set_xticks([-4,-2,0,2,4])
                ax.tick_params(axis='both', labelsize=fs)
                for fmt in fmts: plt.savefig(self.output_path+self.prefix+self.resolution_prefix+'HI_effects_'+model+fmt)

                # write model statistics
                with open(self.output_path+self.prefix+self.resolution_prefix+'parameters_'+model+'.txt', 'w') as ofile:
                    ofile.write('total number of viruses:\t'+str(len(self.viruses))+'\n')
                    ofile.write('number of test viruses:\t'+str(len(self.virus_effect[model]))+'\n')
                    ofile.write('number of reference viruses:\t'+str(len(self.ref_strains))+'\n')
                    ofile.write('number of antisera:\t'+str(len(self.sera))+'\n')
                    ofile.write('number of HI measurements:\t'+str(len(self.HI_normalized))+'\n')
                    ofile.write('number of genetic parameters:\t'+str(self.HI_split_count if model=='tree' else len(self.relevant_muts))+'\n')
                    ofile.write('number of positive genetic parameters:\t'
                                +str(sum([n.dHI>0.001 for n in self.tree.get_nonterminals(order='postorder')]) if model=='tree'
                                     else sum([x>0.001 for x in self.mutation_effects.values()]))+'\n')

            except:
                print ("can't plot effect distributions")


        for map_to_tree, model in [(True, 'tree'), (False,'mutation')]:
            self.map_HI(training_fraction=0.9, method=method,lam_HI=lam_HI, lam_avi=lam_avi,
                        lam_pot = lam_pot, force_redo=True, map_to_tree=map_to_tree, subset_strains=True)

            self.validate(plot=True)
            for fmt in fmts: plt.savefig(self.output_path+self.prefix+self.resolution_prefix+'HI_prediction_virus_'+model+fmt)

            self.map_HI(training_fraction=0.9, method=method,lam_HI=lam_HI, lam_avi=lam_avi,
                        lam_pot = lam_pot, force_redo=True, map_to_tree=map_to_tree)
            self.validate(plot=True)
            for fmt in fmts: plt.savefig(self.output_path+self.prefix+self.resolution_prefix+'HI_prediction_'+model+fmt)

        self.save_trunk_cHI()

    def validate(self, plot=False, cutoff=0.0, validation_set = None, incl_ref_strains='yes'):
        if validation_set is None:
            validation_set=self.test_HI
        from scipy.stats import linregress, pearsonr
        self.validation = {}
        for key, val in validation_set.iteritems():
            if self.map_to_tree:
                pred_HI = self.predict_HI_tree(key[0], key[1], cutoff=cutoff)
            else:
                pred_HI = self.predict_HI_mutations(key[0], key[1], cutoff=cutoff)

            if pred_HI is not None:
                if any([incl_ref_strains=='yes',
                        incl_ref_strains=='no' and (key[0].upper() not in self.ref_strains),
                        incl_ref_strains=='only' and (key[0].upper() in self.ref_strains)]):
                    self.validation[key] = (val, pred_HI)

        a = np.array(self.validation.values())
        print ("number of prediction-measurement pairs",a.shape)
        self.abs_error = np.mean(np.abs(a[:,0]-a[:,1]))
        self.rms_error = np.sqrt(np.mean((a[:,0]-a[:,1])**2))
        self.slope, self.intercept, tmpa, tmpb, tmpc = linregress(a[:,0], a[:,1])
        print ("error (abs/rms): ",self.abs_error, self.rms_error)
        print ("slope, intercept:", self.slope, self.intercept)
        print ("pearson correlation:", pearsonr(a[:,0], a[:,1]))

        if plot:
            import matplotlib.pyplot as plt
            import seaborn as sns
            sns.set_style('darkgrid')
            plt.figure(figsize=(figheight,figheight))
            ax = plt.subplot(111)
            plt.plot([-1,6], [-1,6], 'k')
            plt.scatter(a[:,0], a[:,1])
            plt.ylabel(r"predicted $\log_2$ distance", fontsize = fs)
            plt.xlabel(r"measured $\log_2$ distance" , fontsize = fs)
            ax.tick_params(axis='both', labelsize=fs)
            plt.ylim([-3,8])
            plt.xlim([-3,7])
            plt.text(-2.5,7.3, ('tree model' if self.map_to_tree else 'mutation model'), weight='bold', fontsize=fs)
            plt.text(-2.5,6,'regularization:\nprediction error:', fontsize = fs-2)
            plt.text(1.2,6, str(self.lam_HI)+'/'+str(self.lam_pot)+'/'+str(self.lam_avi)+' (HI/pot/avi)'
                     +'\n'+str(round(self.abs_error, 2))\
                     +'/'+str(round(self.rms_error, 2))+' (abs/rms)', fontsize = fs-2)
            plt.tight_layout()
        return a.shape[0]


    def add_titers(self):
        '''
        write the HI models into the tree structure that is going to be exported
        to auspice for display in the browser
        '''
        for ref in self.ref_strains: # add empty data structures
            self.node_lookup[ref].HI_titers= defaultdict(dict)
            self.node_lookup[ref].HI_titers_raw= defaultdict(dict)
            self.node_lookup[ref].potency_mut={}
            self.node_lookup[ref].potency_tree={}
            self.node_lookup[ref].autologous_titers = {}
        for ref in self.sera: # add serum potencies to reference viruses
            self.node_lookup[ref[0]].autologous_titers[ref[1]] = np.median(self.autologous_titers[ref]['val'])
            if 'mutation' in self.virus_effect:
                self.node_lookup[ref[0]].potency_mut[ref[1]] = self.serum_potency['mutation'][ref]
            if 'tree' in self.virus_effect:
                self.node_lookup[ref[0]].potency_tree[ref[1]] = self.serum_potency['tree'][ref]
        for (test, ref), val in self.HI_normalized.iteritems(): # add normalized HI titers for coloring on tree
            try:
                self.node_lookup[ref[0]].HI_titers[self.node_lookup[test].clade][ref[1]] = val
            except:
                print("Can't assign",test, ref)
        for (test, ref), val in self.HI_raw.iteritems(): # add the raw titer values for display in tool tips
            try:
                self.node_lookup[ref[0]].HI_titers_raw[self.node_lookup[test].clade][ref[1]] = val
            except:
                print("Can't assign",test, ref)
        for test in self.HI_strains: # add the virus avidities to each node with HI information
            if 'mutation' in self.virus_effect:
                self.node_lookup[test].avidity_tree = self.virus_effect['mutation'][test]
            if 'tree' in self.virus_effect:
                self.node_lookup[test].avidity_mut = self.virus_effect['tree'][test]
        for ref in self.ref_strains: # add average potencies and titers to each reference virus (not actually used in auspice right now).
            self.node_lookup[ref].mean_HI_titers = {key:np.mean(titers.values()) for key, titers in
                                                self.node_lookup[ref].HI_titers.iteritems()}
            self.node_lookup[ref].mean_potency_tree = np.mean(self.node_lookup[ref].potency_tree.values())
            self.node_lookup[ref].mean_potency_mut = np.mean(self.node_lookup[ref].potency_mut.values())

    def check_sources(self):
        self.source_HIs = defaultdict(dict)
        with myopen(self.HI_fname, 'r') as infile:
            for line in infile:
                test, ref_virus, serum, src_id, val_str = line.strip().split()
                try:
                    val = float(val_str)
                    if not np.isnan(val):
                        self.source_HIs[src_id][test, (ref_virus, serum)] = self.normalize((ref_virus, serum), float(val))
                except:
                    print (test, ref_virus, serum, src_id, float(val))

        self.source_validation = {}
        for src_id in self.source_HIs:
            print ('\n############### \n',src_id,'\n############### \n')
            print ("number of measurements:",len(self.source_HIs[src_id]))
            try:
                n_checks = self.validate(validation_set=self.source_HIs[src_id], incl_ref_strains='no')
                self.source_validation[src_id] = [self.abs_error, self.rms_error, self.slope, self.intercept, n_checks]
            except:
                print ("skipped due to too few measurements")

    def save_trunk_cHI(self):
        cHI_trunk = []
        trunk_muts = []
        co=0.5
        tmp_pivots = self.tree.root.pivots
        for node in self.tree.get_nonterminals(order='postorder'):
            node.num_date = np.min([c.num_date for c in node.clades])
            if node.trunk:
                tmp_freq = node.freq['global']
                if tmp_freq[0]>0.5:
                    continue
                ii = np.argmax(tmp_freq>co)
                slope = (tmp_freq[ii]-tmp_freq[ii-1])/(tmp_pivots[ii]-tmp_pivots[ii-1])
                tp = tmp_pivots[ii-1] + (co-tmp_freq[ii-1])/slope

                tmp_muts = [x for x in node.aa_muts.items() if x[1]!='']
                for muts in tmp_muts:
                    gene = muts[0]
                    for pos in muts[1].split(','):
                        tmp_mut = (gene, pos)
                        if tmp_mut in self.mutation_effects:
                            trunk_muts.append([tp,self.mutation_effects[tmp_mut]])
                cHI_trunk.append([tp, node.cHI] )

        cHI_trunk = np.array(cHI_trunk)[::-1]
        trunk_muts = np.array(trunk_muts)[::-1]
        np.savetxt('data/'+self.prefix+self.resolution+'_cHI.txt', cHI_trunk)
        np.savetxt('data/'+self.prefix+self.resolution+'_trunk_muts.txt', trunk_muts)

        from random import sample
        n_leafs = 1
        leaf_sample = []
        for y in range(int(self.time_interval[0]), int(self.time_interval[1])):
            tmp_leafs = [leaf for leaf in self.tree.get_terminals() if leaf.num_date>y and leaf.num_date<y+1]
            leaf_sample.extend(sample(tmp_leafs, min(n_leafs, len(tmp_leafs))))

        cHI_trunk = []
        for node in leaf_sample:
            tmp = []
            while node.parent_node is not None:
                tmp.append([node.num_date, node.cHI])
                node = node.parent_node
            cHI_trunk.append(np.array(tmp)[::-1])
        import cPickle as pickle
        with open('data/'+self.prefix+self.resolution+'_cHI_path.pkl', 'w') as ofile:
            pickle.dump(cHI_trunk, ofile)

    def predict_HI_tree(self, virus, serum, cutoff=0.0):
        path = self.get_path_no_terminals(virus,serum[0])
        if path is not None:
            return self.serum_potency['tree'][serum] \
                    + self.virus_effect['tree'][virus] \
                    + np.sum([b.dHI for b in path if b.dHI>cutoff])
        else:
            return None

    def predict_HI_mutations(self, virus, serum, cutoff=0.0):
        muts= self.get_mutations(serum[0], virus)
        if muts is not None:
            return self.serum_potency['mutation'][serum] \
                    + self.virus_effect['mutation'][virus] \
                    + np.sum([self.mutation_effects[mut] for mut in muts
                    if (mut in self.mutation_effects and self.mutation_effects[mut]>cutoff)])
        else:
            return None

