import os
import logging
import numpy as np
import time
from collections import defaultdict
import pandas as pd
from pprint import pprint
import sys

from .utils import myopen

TITER_ROUND=4

class InsufficientDataException(Exception):
    pass

class TiterCollection(object):
    """
    Container for raw titer values and methods for analyzing these values.
    """
    @staticmethod
    def load_from_file(filenames, excluded_sources=None):
        """Load titers from a tab-delimited file.

        Parameters
        ----------
        filename : str
            tab-delimited file containing titer strains, serum,
            and values
        excluded_sources : list of str
            sources in the titers file to exclude

        Returns
        -------
        tuple (dict, list, list)
            tuple of a dict of titer measurements, list of strains, list of sources

        >>> measurements, strains, sources = TiterCollection.load_from_file("tests/data/titer_model/h3n2_titers_subset.tsv")
        >>> type(measurements)
        <class 'dict'>
        >>> len(measurements)
        11
        >>> measurements[("A/Acores/11/2013", ("A/Alabama/5/2010", "F27/10"))]
        [80.0]
        >>> len(strains)
        13
        >>> len(sources)
        5
        >>> measurements, strains, sources = TiterCollection.load_from_file("tests/data/titer_model/h3n2_titers_subset.tsv", excluded_sources=["NIMR_Sep2013_7-11.csv"])
        >>> len(measurements)
        5
        >>> measurements.get(("A/Acores/11/2013", ("A/Alabama/5/2010", "F27/10")))
        >>>
        >>> output = TiterCollection.load_from_file("tests/data/titer_model/missing.tsv")
        Traceback (most recent call last):
          File "<ipython-input-2-0ea96a90d45d>", line 1, in <module>
            open("tests/data/titer_model/missing.tsv", "r")
        FileNotFoundError: [Errno 2] No such file or directory: 'tests/data/titer_model/missing.tsv'
        """
        if excluded_sources is None:
            excluded_sources = []

        measurements = defaultdict(list)
        strains = set()
        sources = set()
        titer_files = [filenames] if type(filenames)==str else filenames

        for fname in titer_files:
            with myopen(fname, 'r') as infile:
                for line in infile:
                    entries = line.strip().split('\t')
                    try:
                        val = float(entries[4])
                    except:
                        continue
                    test, ref_virus, serum, src_id = (entries[0], entries[1],entries[2],
                                                      entries[3])

                    ref = (ref_virus, serum)
                    if src_id not in excluded_sources:
                        try:
                            measurements[(test, (ref_virus, serum))].append(val)
                            strains.update([test, ref_virus])
                            sources.add(src_id)
                        except:
                            print(line.strip())

        print("Read titers from %s, found:" % ' '.join(titer_files), file=sys.stderr)
        print(" --- %i strains" % len(strains), file=sys.stderr)
        print(" --- %i data sources" % len(sources), file=sys.stderr)
        print(" --- %i total measurements" % sum([len(x) for x in measurements.values()]), file=sys.stderr)

        return dict(measurements), list(strains), list(sources)

    @staticmethod
    def count_strains(titers):
        """Count test and reference virus strains in the given titers.

        Parameters
        ----------
        titers : defaultdict
            titer measurements indexed by test, reference,
            and serum

        Returns
        -------
        dict
            number of measurements per strain

        >>> measurements, strains, sources = TiterCollection.load_from_file("tests/data/titer_model/h3n2_titers_subset.tsv")
        >>> titer_counts = TiterCollection.count_strains(measurements)
        >>> titer_counts["A/Acores/11/2013"]
        6
        >>> titer_counts["A/Acores/SU43/2012"]
        3
        >>> titer_counts["A/Cairo/63/2012"]
        2
        """
        counts = defaultdict(int)
        for key in titers.keys():
            measurements = len(titers[key])
            counts[key[0]] += measurements

        return counts

    @staticmethod
    def filter_strains(titers, strains):
        """Filter the given titers to only include values from the given strains
        (test or reference).

        Parameters
        ----------
        titers : dict
            titer values indexed by test and reference strain and
            serum
        strains : list
            names of strains to keep titers for

        Returns
        -------
        dict
            reduced dictionary of titer measurements containing only those were
            test and reference virus are part of the strain list

        >>> measurements, strains, sources = TiterCollection.load_from_file("tests/data/titer_model/h3n2_titers_subset.tsv")
        >>> len(measurements)
        11

        Test the case when a test strain exists in the subset but the none of
        its corresponding reference strains do.

        >>> len(TiterCollection.filter_strains(measurements, ["A/Acores/11/2013"]))
        0

        Test when both the test and reference strains exist in the subset.

        >>> len(TiterCollection.filter_strains(measurements, ["A/Acores/11/2013", "A/Alabama/5/2010", "A/Athens/112/2012"]))
        2
        >>> len(TiterCollection.filter_strains(measurements, ["A/Acores/11/2013", "A/Acores/SU43/2012", "A/Alabama/5/2010", "A/Athens/112/2012"]))
        3
        >>> len(TiterCollection.filter_strains(measurements, []))
        0
        """
        return {key: value for key, value in titers.items()
                if key[0] in strains and key[1][0] in strains}


    def __init__(self, titers, **kwargs):
        """Accepts the name of a file containing titers to load or a preloaded titers
        dictionary.

        Parameters
        ----------
        titers : TYPE
            Description
        **kwargs
            Description
        """
        self.kwargs = kwargs

        # Assign titers and prepare list of strains.
        if (isinstance(titers, str) and os.path.isfile(titers))\
            or isinstance(titers, list):
            self.read_titers(titers)
        else:
            self.titers = titers
            strain_counts = type(self).count_strains(titers)
            self.strains = strain_counts.keys()

    def read_titers(self, fname):
        self.titer_fname = fname
        if "excluded_tables" in self.kwargs:
            self.excluded_tables = self.kwargs["excluded_tables"]
        else:
            self.excluded_tables = []

        self.titers, self.strains, self.sources = type(self).load_from_file(
                                                            fname,
                                                            self.excluded_tables
                                                        )

    def normalize(self, ref, val):
        '''
        take the log2 difference of test titers and autologous titers

        Parameters
        ----------
        ref : TYPE
            Description
        val : TYPE
            Description

        Returns
        -------
        TYPE
            Description
        '''
        consensus_func = np.mean
        return consensus_func(np.log2(self.autologous_titers[ref]['val'])) \
                - consensus_func(np.log2(val))

    def determine_autologous_titers(self):
        '''
        scan the titer measurements for autologous (self) titers and make a dictionary
        stored in self to look them up later. If no autologous titer is found, use the
        maximum titer. This follows the rationale that test titers are generally lower
        than autologous titers and the highest test titer is often a reasonably
        approximation of the autologous titer.
        '''
        autologous = defaultdict(list)
        all_titers_per_serum = defaultdict(list)
        for (test, ref), val in self.titers.items():
            all_titers_per_serum[ref].append(val)
            if ref[0]==test:
                autologous[ref].append(val)

        self.autologous_titers = {}
        for serum in all_titers_per_serum:
            if serum in autologous:
                self.autologous_titers[serum] = {'val':autologous[serum], 'autologous':True}
                #print("autologous titer found for",serum)
            else:
                # use max tier if there are at least 10 measurements, don't bother otherwuise
                if len(all_titers_per_serum[serum])>10:
                    autologous_proxy = np.percentile([np.median(x) for x in all_titers_per_serum[serum]],90)
                    self.autologous_titers[serum] = {'val':autologous_proxy,
                                                     'autologous':False}
                    print(serum,": using 90% percentile instead of autologous,",
                          autologous_proxy)
                else:
                    pass
                    # print("discarding",serum,"since there are only ",
                    #        len(all_titers_per_serum[serum]),'measurements')

    def normalize_titers(self):
        '''
        convert the titer measurements into the log2 difference between the average
        titer measured between test virus and reference serum and the average
        homologous titer. all measurements relative to sera without homologous titer
        are excluded
        '''
        self.determine_autologous_titers()

        self.titers_normalized = {}
        self.consensus_titers_raw = {}
        self.measurements_per_serum = defaultdict(int)
        for (test, ref), val in self.titers.items():
            if ref in self.autologous_titers: # use only titers for which estimates of the autologous titer exists
                self.titers_normalized[(test, ref)] = self.normalize(ref, val)
                self.consensus_titers_raw[(test, ref)] = np.median(val)
                self.measurements_per_serum[ref]+=1
            else:
                pass
                #print("no homologous titer found:", ref)

    def strain_census(self, titers):
        """
        make lists of reference viruses, test viruses and sera
        (there are often multiple sera per reference virus)

        >>> measurements, strains, sources = TiterCollection.load_from_file("tests/data/titer_model/h3n2_titers_subset.tsv")
        >>> titers = TiterCollection(measurements)
        >>> sera, ref_strains, test_strains = titers.strain_census(measurements)
        >>> len(sera)
        9
        >>> len(ref_strains)
        9
        >>> len(test_strains)
        13

        Parameters
        ----------
        titers : TYPE
            Description

        Returns
        -------
        TYPE
            Description
        """
        sera = set()
        ref_strains = set()
        test_strains = set()

        for test, ref in titers:
            test_strains.add(test)
            test_strains.add(ref[0])
            sera.add(ref)
            ref_strains.add(ref[0])

        return list(sera), list(ref_strains), list(test_strains)


class TiterModel(object):
    '''
    this class fits a linear model to titer measurements using
    different models that describe titer differences in a parsimonious way.
    Two additive models are currently implemented, the tree and the substitution
    model. The tree model describes titer drops as a sum of terms associated with
    branches in the tree, while the substitution model attributes titer drops to amino
    acid mutations. More details on the methods can be found in
    Neher et al, PNAS, 2016

    '''
    def __init__(self, serum_Kc=0, **kwargs):
        '''
        default constructor assumes a Bio.Phylo tree as first positional argument and a dictionay of
        titer measurements as second positional argument. This dictionary has composite keys consisting
        of the (test_virus_strain_name, (reference_virus_strain_name, serum_id))

        Parameters
        ----------
        serum_Kc : int, optional
            optional argument that can be used to even out contribution of sera.
            should be roughly the inverse of the number of measurements beyond which
            the contribution of a serum should saturate
        **kwargs
            other keyword arguments
        '''
        self.kwargs = kwargs
        self.serum_Kc = serum_Kc


    def assign_titers(self, titers, strains):
        # Load titer measurements from a file or from a given dictionary of
        # measurements.
        if (isinstance(titers, str) and os.path.isfile(titers))\
            or isinstance(titers, list):
            titer_measurements, in_strains, sera = TiterCollection.load_from_file(titers)
        else:
            titer_measurements = titers

        # Filter titer measurements to those from strains in the strain list.
        filtered_titer_measurements = TiterCollection.filter_strains(
            titer_measurements,
            strains
        )

        # Create a titer collection for the filtered titer measurements.
        self.titers = TiterCollection(filtered_titer_measurements)

        # Normalize titers.
        self.titers.normalize_titers()

        # Determine distinct sera, reference strains, and test strains.
        self.sera, self.ref_strains, self.test_strains = self.titers.strain_census(self.titers.titers_normalized)
        print("Normalized titers and restricted to measurements in list:")
        self.titer_stats()


    def titer_stats(self):
        print(" - remaining data set")
        print(' ---', len(self.ref_strains), " reference virues")
        print(' ---', len(self.sera), " sera")
        print(' ---', len(self.test_strains), " test_viruses")
        print(' ---', len(self.titers.titers_normalized), " non-redundant test virus/serum pairs")
        if hasattr(self, 'train_titers'):
            print(' ---', len(self.train_titers), " measurements in training set")


    def make_training_set(self, training_fraction=1.0, subset_strains=False, **kwargs):
        if training_fraction<1.0: # validation mode, set aside a fraction of measurements to validate the fit
            self.test_titers, self.train_titers = {}, {}
            if subset_strains:    # exclude a fraction of test viruses as opposed to a fraction of the titers
                from random import sample
                tmp = set(self.test_strains)
                tmp.difference_update(self.ref_strains) # don't use references viruses in the set to sample from
                training_strains = sample(tmp, int(training_fraction*len(tmp)))
                for tmpstrain in self.ref_strains:      # add all reference viruses to the training set
                    if tmpstrain not in training_strains:
                        training_strains.append(tmpstrain)
                for key, val in self.titers.titers_normalized.items():
                    if key[0] in training_strains:
                        self.train_titers[key]=val
                    else:
                        self.test_titers[key]=val
            else: # simply use a fraction of all measurements for testing
                for key, val in self.titers.titers_normalized.items():
                    if np.random.uniform()>training_fraction:
                        self.test_titers[key]=val
                    else:
                        self.train_titers[key]=val
        else: # without the need for a test data set, use the entire data set for training
            self.train_titers = self.titers.titers_normalized

        self.sera, self.ref_strains, self.test_strains = self.titers.strain_census(self.train_titers)
        print("Made training data as fraction",training_fraction, "of all measurements")
        self.titer_stats()


    def _train(self, method='nnl1reg',  lam_drop=1.0, lam_pot = 0.5, lam_avi = 3.0, **kwargs):
        '''
        determine the model parameters -- lam_drop, lam_pot, lam_avi are
        the regularization parameters.

        Parameters
        ----------
        method : str, optional
            Description
        lam_drop : float, optional
            Description
        lam_pot : float, optional
            Description
        lam_avi : float, optional
            Description
        **kwargs
            Description
        '''
        self.lam_pot = lam_pot
        self.lam_avi = lam_avi
        self.lam_drop = lam_drop
        if len(self.train_titers)==0:
            raise InsufficientDataException("Error: No titers in training set.")
        else:
            if method=='l1reg':  # l1 regularized fit, no constraint on sign of effect
                self.model_params = self.fit_l1reg()
            elif method=='nnls':  # non-negative least square, not regularized
                self.model_params = self.fit_nnls()
            elif method=='nnl2reg': # non-negative L2 norm regularized fit
                self.model_params = self.fit_nnl2reg()
            elif method=='nnl1reg':  # non-negative fit, branch terms L1 regularized, avidity terms L2 regularized
                self.model_params = self.fit_nnl1reg()

            print('rms deviation on training set=',np.sqrt(self.fit_func()))

        # extract and save the potencies and virus effects. The genetic parameters
        # are subclass specific and need to be process by the subclass
        self.serum_potency = {serum:self.model_params[self.genetic_params+ii]
                              for ii, serum in enumerate(self.sera)}
        self.virus_effect = {strain:self.model_params[self.genetic_params+len(self.sera)+ii]
                             for ii, strain in enumerate(self.test_strains)}


    def fit_func(self):
        return np.mean( (self.titer_dist - np.dot(self.design_matrix, self.model_params))**2 )


    def validate(self, plot=False, cutoff=0.0, validation_set = None, fname=None):
        '''
        predict titers of the validation set (separate set of test_titers aside previously)
        and compare against known values. If requested by plot=True,
        a figure comparing predicted and measured titers is produced

        Compute basic error metrics for actual vs. predicted titer values.
        Return a dictionary of {'metric': computed_metric, 'values': [(actual, predicted), ...]}, save a copy in self.validation

        Parameters
        ----------
        plot : bool, optional
            Description
        cutoff : float, optional
            Description
        validation_set : None, optional
            Description
        fname : None, optional
            Description

        Returns
        -------
        TYPE
            Description
        '''
        from scipy.stats import linregress, pearsonr
        if validation_set is None:
            validation_set=self.test_titers
        validation = {}
        for key, val in validation_set.items():
            pred_titer = self.predict_titer(key[0], key[1], cutoff=cutoff)
            validation[key] = (val, pred_titer)

        validation_array = np.array(validation.values())
        actual = validation_array[:,0]
        predicted = validation_array[:,1]

        regression = linregress(actual, predicted)
        model_performance = {
                        'slope': regression[0],
                        'intercept': regression[1],
                        'r_squared': pearsonr(actual, predicted)[0]**2,
                        'abs_error':  np.mean(np.abs(actual-predicted)),
                        'rms_error': np.sqrt(np.mean((actual-predicted)**2)),
        }
        pprint(model_performance)
        model_performance['values'] = validation.values()

        self.validation = model_performance

        if plot:
            try:
                import matplotlib.pyplot as plt
                import seaborn as sns
            except ImportError:
                print("Plotting requires a working matplotlib and seaborn installation.")
            else:
                fs=16
                sns.set_style('darkgrid')
                plt.figure()
                ax = plt.subplot(111)
                plt.plot([-1,6], [-1,6], 'k')
                plt.scatter(actual, predicted)
                plt.ylabel(r"predicted $\log_2$ distance", fontsize = fs)
                plt.xlabel(r"measured $\log_2$ distance" , fontsize = fs)
                ax.tick_params(axis='both', labelsize=fs)
                plt.text(-2.5,6,'regularization:\nprediction error:\nR^2:', fontsize = fs-2)
                plt.text(1.2,6, str(self.lam_drop)+'/'+str(self.lam_pot)+'/'+str(self.lam_avi)+' (HI/pot/avi)'
                         +'\n'+str(round(model_performance['abs_error'], 2))+'/'+str(round(model_performance['rms_error'], 2))+' (abs/rms)'
                         + '\n' + str(model_performance['r_squared']), fontsize = fs-2)
                plt.tight_layout()

                if fname:
                    plt.savefig(fname)

        return model_performance

    def reference_virus_statistic(self):
        '''
        count measurements for every reference virus and serum
        '''
        def dstruct():
            return defaultdict(int)
        self.titer_counts = defaultdict(dstruct)
        for test_vir, (ref_vir, serum) in self.titers.titers_normalized:
            self.titer_counts[ref_vir][serum]+=1


    def compile_titers(self):
        '''
        compiles titer measurements into a json file organized by reference virus
        during visualization, we need the average distance of a test virus from
        a reference virus across sera. hence the hierarchy [ref][test][serum]
        NOTE: this uses node.name instead of node.clade

        Returns
        -------
        TYPE
            Description
        '''
        def dstruct():
            return defaultdict(dict)
        titer_json = defaultdict(dstruct)

        for key, val in self.titers.titers_normalized.items():
            test_clade, (ref_clade, serum) = key
            titer_json[ref_clade][test_clade][serum] = [np.round(val,TITER_ROUND), np.median(self.titers.titers[key])]

        return titer_json


    def compile_potencies(self):
        '''
        compile a json structure containing potencies for visualization
        we need rapid access to all sera for a given reference virus, hence
        the structure is organized by [ref][serum]

        Returns
        -------
        TYPE
            Description
        '''
        potency_json = defaultdict(dict)
        for (ref_clade, serum), val in self.serum_potency.items():
            potency_json[ref_clade][serum] = np.round(val,TITER_ROUND)

        # add the average potency (weighed by the number of measurements per serum)
        # to the exported data structure
        self.reference_virus_statistic()
        mean_potency = defaultdict(int)
        for (ref_vir, serum), val in self.serum_potency.items():
            mean_potency[ref_vir] += self.titer_counts[ref_vir][serum]*val
        for ref_vir in self.ref_strains:
            potency_json[ref_vir]['mean_potency'] = 1.0*mean_potency[ref_vir]/np.sum(list(self.titer_counts[ref_vir].values()))

        return potency_json


    def compile_virus_effects(self):
        '''
        compile a json structure containing virus_effects for visualization

        Returns
        -------
        TYPE
            Description
        '''
        return {test_vir:np.round(val,TITER_ROUND) for test_vir, val in self.virus_effect.items()}


    ##########################################################################################
    # define fitting routines for different objective functions
    ##########################################################################################
    def fit_l1reg(self):
        '''
        regularize genetic parameters with an l1 norm regardless of sign

        Returns
        -------
        TYPE
            Description
        '''
        try:
            from cvxopt import matrix, solvers
        except ImportError:
            raise ImportError("To infer titer models, you need a working installation of cvxopt")
        n_params = self.design_matrix.shape[1]
        n_genetic = self.genetic_params
        n_sera = len(self.sera)
        n_v = len(self.test_strains)

        # set up the quadratic matrix containing the deviation term (linear xterm below)
        # and the l2-regulatization of the avidities and potencies
        P1 = np.zeros((n_params+n_genetic,n_params+n_genetic))
        P1[:n_params, :n_params] = self.TgT
        for ii in range(n_genetic, n_genetic+n_sera):
            P1[ii,ii]+=self.lam_pot
        for ii in range(n_genetic+n_sera, n_params):
            P1[ii,ii]+=self.lam_avi
        P = matrix(P1)

        # set up cost for auxillary parameter and the linear cross-term
        q1 = np.zeros(n_params+n_genetic)
        q1[:n_params] = -np.dot( self.titer_dist, self.design_matrix)
        q1[n_params:] = self.lam_drop
        q = matrix(q1)

        # set up linear constraint matrix to regularize the HI parametesr
        h = matrix(np.zeros(2*n_genetic))   # Gw <=h
        G1 = np.zeros((2*n_genetic,n_params+n_genetic))
        G1[:n_genetic, :n_genetic] = -np.eye(n_genetic)
        G1[:n_genetic:, n_params:] = -np.eye(n_genetic)
        G1[n_genetic:, :n_genetic] = np.eye(n_genetic)
        G1[n_genetic:, n_params:] = -np.eye(n_genetic)
        G = matrix(G1)
        W = solvers.qp(P,q,G,h)
        return np.array([x for x in W['x']])[:n_params]


    def fit_nnls(self):
        from scipy.optimize import nnls
        return nnls(self.design_matrix, self.titer_dist)[0]


    def fit_nnl2reg(self):
        try:
            from cvxopt import matrix, solvers
        except ImportError:
            raise ImportError("To infer titer models, you need a working installation of cvxopt")
        n_params = self.design_matrix.shape[1]
        P = matrix(np.dot(self.design_matrix.T, self.design_matrix) + self.lam_drop*np.eye(n_params))
        q = matrix( -np.dot( self.titer_dist, self.design_matrix))
        h = matrix(np.zeros(n_params)) # Gw <=h
        G = matrix(-np.eye(n_params))
        W = solvers.qp(P,q,G,h)
        return np.array([x for x in W['x']])


    def fit_nnl1reg(self):
        '''l1 regularization of titer drops with non-negativity constraints

        Returns
        -------
        TYPE
            Description
        '''
        try:
            from cvxopt import matrix, solvers
        except ImportError:
            raise ImportError("To infer titer models, you need a working installation of cvxopt")
        n_params = self.design_matrix.shape[1]
        n_genetic = self.genetic_params
        n_sera = len(self.sera)
        n_v = len(self.test_strains)

        # set up the quadratic matrix containing the deviation term (linear xterm below)
        # and the l2-regulatization of the avidities and potencies
        P1 = np.zeros((n_params,n_params))
        P1[:n_params, :n_params] = self.TgT
        for ii in range(n_genetic, n_genetic+n_sera):
            P1[ii,ii]+=self.lam_pot
        for ii in range(n_genetic+n_sera, n_params):
            P1[ii,ii]+=self.lam_avi
        P = matrix(P1)

        # set up cost for auxillary parameter and the linear cross-term
        q1 = np.zeros(n_params)
        q1[:n_params] = -np.dot(self.titer_dist, self.design_matrix)
        q1[:n_genetic] += self.lam_drop
        q = matrix(q1)

        # set up linear constraint matrix to enforce positivity of the
        # dTiters and bounding of dTiter by the auxillary parameter
        h = matrix(np.zeros(n_genetic))     # Gw <=h
        G1 = np.zeros((n_genetic,n_params))
        G1[:n_genetic, :n_genetic] = -np.eye(n_genetic)
        G = matrix(G1)
        W = solvers.qp(P,q,G,h)
        return np.array([x for x in W['x']])[:n_params]

##########################################################################################
# END GENERIC CLASS
##########################################################################################



##########################################################################################
# TREE MODEL
##########################################################################################
class TreeModel(TiterModel):
    """
    tree_model extends titers and fits the antigenic differences
    in terms of contributions on the branches of the phylogenetic tree.
    nodes in the tree are decorated with attributes 'dTiter' that contain
    the estimated titer drops across the branch
    """
    def __init__(self, tree, titers, *args, **kwargs):
        super(TreeModel, self).__init__(*args, **kwargs)
        self.tree = None
        self.prepare_tree(tree)
        strains = [x.name for x in self.tree.get_terminals()]
        self.assign_titers(titers, strains)


    def prepare_tree(self, tree):
        self.tree = tree # not copied, just linked
        # produce dictionaries that map node names to nodes regardless of capitalization
        self.strain_lookup = {n.name:n for n in tree.get_terminals()}
        self.strain_lookup.update({n.name.upper():n for n in tree.get_terminals()})
        self.strain_lookup.update({n.name.lower():n for n in tree.get_terminals()})

        # have each node link to its parent. this will be needed for walking up and down the tree
        # but should be already in place if treetime is used.
        self.tree.root.up=None
        for node in self.tree.get_nonterminals():
            for c in node.clades:
                c.up = node



    def cross_validate(self, n, **kwargs):
        '''
        For each of n iterations, randomly re-allocate titers to training and test set.
        Fit the model using training titers, assess performance using test titers (see TiterModel.validate)
        Append dictionaries of {'abs_error': , 'rms_error': , 'values': [(actual, predicted), ...], etc.} for each iteration to the model_performance list.
        Return model_performance, and save a copy in self.cross_validation

        Parameters
        ----------
        n : TYPE
            Description
        **kwargs
            Description

        Returns
        -------
        TYPE
            Description
        '''

        model_performance = []
        for iteration in range(n):
            self.prepare(**kwargs) # randomly reassign titers to training and test sets
            self.train(**kwargs) # train the model
            performance = self.validate() # assess performance on the withheld test data. Returns {'values': [(actual, predicted), ...], 'metric': metric_value, ...}
            model_performance.append(performance)

        self.cross_validation = model_performance
        return self.cross_validation

    def prepare(self, **kwargs):
        self.make_training_set(**kwargs)
        self.find_titer_splits(criterium= kwargs['criterium']
                                          if 'criterium' in kwargs else None)
        if len(self.train_titers)>1:
            self.make_treegraph()
        else:
            raise InsufficientDataException("TreeModel: Not enough titers in training set, found {}".format(len(self.train_titers)))

    def get_path_no_terminals(self, v1, v2):
        '''
        returns the path between two tips in the tree excluding the terminal branches.

        Parameters
        ----------
        v1 : TYPE
            Description
        v2 : TYPE
            Description

        Returns
        -------
        TYPE
            Description
        '''
        if v1 in self.strain_lookup and v2 in self.strain_lookup:
            p1 = [self.strain_lookup[v1]]
            p2 = [self.strain_lookup[v2]]
            for tmp_p in [p1,p2]:
                while tmp_p[-1].up != self.tree.root:
                    tmp_p.append(tmp_p[-1].up)
                tmp_p.append(self.tree.root)
                tmp_p.reverse()

            for pi, (tmp_v1, tmp_v2) in enumerate(zip(p1,p2)):
                if tmp_v1!=tmp_v2:
                    break
            path = [n for n in p1[pi:] if n.titer_info>1] + [n for n in p2[pi:] if n.titer_info>1]
        else:
            path = None
        return path


    def find_titer_splits(self, criterium=None):
        '''
        walk through the tree, mark all branches that are to be included as model variables
         - no terminals
         - criterium: callable that can be used to exclude branches e.g. if
                      amino acid mutations map to this branch.

        Parameters
        ----------
        criterium : None, optional
            Description
        '''
        if criterium is None:
            criterium = lambda x:True
        # flag all branches on the tree with titer_info = True if they lead to strain with titer data
        for leaf in self.tree.get_terminals():
            if leaf.name in self.test_strains:
                leaf.serum = leaf.name in self.ref_strains
                leaf.titer_info = 1
            else:
                leaf.serum, leaf.titer_info=False, 0

        for node in self.tree.get_nonterminals(order='postorder'):
            node.titer_info = sum([c.titer_info for c in node.clades])
            node.serum= False

        # combine sets of branches that span identical sets of titers
        self.titer_split_count = 0  # titer split counter
        self.titer_split_to_branch = defaultdict(list)
        for node in self.tree.find_clades(order='preorder'):
            node.dTiter, node.cTiter, node.constraints = 0, 0, 0
            if node.titer_info>1 and criterium(node):
                node.titer_branch_index = self.titer_split_count
                self.titer_split_to_branch[node.titer_branch_index].append(node)
                # at a bi- or multifurcation, increase the split count and HI index
                # either individual child branches have enough HI info be counted,
                # or the pre-order node iteraction will move towards the root
                if sum([c.titer_info>0 for c in node.clades])>1:
                    self.titer_split_count+=1
                elif node.is_terminal():
                    self.titer_split_count+=1
            else:
                node.titer_branch_index=None

        self.genetic_params = self.titer_split_count
        print("# of reference strains:",len(self.sera))
        print("# of eligible branches with titer constraints", self.titer_split_count)


    def make_treegraph(self):
        '''
        code the path between serum and test virus of each HI measurement into a matrix
        the matrix has dimensions #measurements x #tree branches with HI info
        if the path between test and serum goes through a branch,
        the corresponding matrix element is 1, 0 otherwise
        '''
        tree_graph = []
        titer_dist = []
        weights = []
        # mark HI splits have to have been run before, assigning self.titer_split_count
        n_params = self.titer_split_count + len(self.sera) + len(self.test_strains)
        for (test, ref), val in self.train_titers.items():
            if not np.isnan(val) and not np.isinf(val):
                try:
                    if ref[0] in self.strain_lookup and test in self.strain_lookup:
                        path = self.get_path_no_terminals(test, ref[0])
                        tmp = np.zeros(n_params, dtype=int)
                        # determine branch indices on path
                        branches = np.unique([c.titer_branch_index for c in path
                                                 if c.titer_branch_index is not None])

                        if len(branches): tmp[branches] = 1
                        # add serum effect for heterologous viruses
                        if ref[0]!=test:
                            tmp[self.titer_split_count+self.sera.index(ref)] = 1
                        # add virus effect
                        tmp[self.titer_split_count+len(self.sera)+self.test_strains.index(test)] = 1
                        # append model and fit value to lists tree_graph and titer_dist
                        tree_graph.append(tmp)
                        titer_dist.append(val)
                        weights.append(1.0/(1.0 + self.serum_Kc*self.titers.measurements_per_serum[ref]))
                except:
                    import ipdb; ipdb.set_trace()
                    print(test, ref, "ERROR")

        # convert to numpy arrays and save product of tree graph with its transpose for future use
        self.weights = np.sqrt(weights)
        self.titer_dist =  np.array(titer_dist)*self.weights
        self.design_matrix = (np.array(tree_graph).T*self.weights).T
        self.TgT = np.dot(self.design_matrix.T, self.design_matrix)
        print("Found", self.design_matrix.shape, "measurements x parameters")

    def train(self,**kwargs):
        self._train(**kwargs)
        for node in self.tree.find_clades(order='postorder'):
            node.dTiter=0 # reset branch properties -- only neede for tree model
            node.cTiter=0
        for titer_split, branches in self.titer_split_to_branch.items():
            likely_branch = branches[np.argmax([b.branch_length for b in branches])]
            likely_branch.dTiter = self.model_params[titer_split]
            likely_branch.constraints = self.design_matrix[:,titer_split].sum()

        # integrate the tree model dTiter into a cumulative antigentic evolution score cTiter
        for node in self.tree.find_clades(order='preorder'):
            if node.up is not None:
                node.cTiter = node.up.cTiter + node.dTiter
            else:
                node.cTiter=0

    def predict_titer(self, virus, serum, cutoff=0.0):
        path = self.get_path_no_terminals(virus,serum[0])
        if path is not None:
            pot = self.serum_potency[serum] if serum in self.serum_potency else 0.0
            avi = self.virus_effect[virus] if virus in self.virus_effect else 0.0
            return avi + pot + np.sum([b.dTiter for b in path if b.dTiter>cutoff])
        else:
            return None



##########################################################################################
# SUBSTITUTION MODEL
##########################################################################################
class SubstitutionModel(TiterModel):
    """
    substitution_model extends titers and implements a model that
    seeks to describe titer differences by sums of contributions of
    substitions separating the test and reference viruses. Sequences
    are assumed to be attached to each terminal node in the tree as
    node.translations
    """
    def __init__(self, alignments, titers, *args, **kwargs):
        super(SubstitutionModel, self).__init__(*args, **kwargs)
        self.sequences = defaultdict(dict)
        self.proteins = list(alignments.keys())
        self.substitution_effect={}
        for gene, aln in alignments.items():
            for x in aln:
                self.sequences[x.name][gene] = str(x.seq)

        strains = list(self.sequences.keys())
        self.assign_titers(titers, strains)


    def prepare(self, **kwargs):
        self.make_training_set(**kwargs)
        self.determine_relevant_mutations()
        if len(self.train_titers)>1:
            self.make_seqgraph()
        else:
            raise InsufficientDataException("SubstitutionModel: Not enough titers in training set, found {}".format(len(self.train_titers)))


    def get_mutations(self, strain1, strain2):
        '''return amino acid mutations between viruses specified by strain names as tuples (HA1, F159S)

        Parameters
        ----------
        strain1 : TYPE
            Description
        strain2 : TYPE
            Description

        Returns
        -------
        TYPE
            Description
        '''
        if strain1 in self.sequences and strain2 in self.sequences:
            muts = []
            for prot in self.proteins:
                seq1 = self.sequences[strain1][prot]
                seq2 = self.sequences[strain2][prot]
                muts.extend([(prot, aa1+str(pos+1)+aa2) for pos, (aa1, aa2)
                            in enumerate(zip(seq1, seq2)) if aa1!=aa2])
            return muts
        else:
            return None


    def determine_relevant_mutations(self, min_count=10):
        # count how often each mutation separates a reference test virus pair
        self.mutation_counter = defaultdict(int)
        for (test, ref), val in self.train_titers.items():
            muts = self.get_mutations(ref[0], test)
            if muts is None:
                continue
            for mut in muts:
                self.mutation_counter[mut]+=1

        # make a list of mutations deemed relevant via frequency thresholds
        relevant_muts = []
        for mut, count in self.mutation_counter.items():
            gene = mut[0]
            pos = int(mut[1][1:-1])-1
            aa1, aa2 = mut[1][0],mut[1][-1]
            if count>min_count:
                relevant_muts.append(mut)

        relevant_muts.sort() # sort by gene
        relevant_muts.sort(key = lambda x:int(x[1][1:-1])) # sort by position in gene
        self.relevant_muts = relevant_muts
        self.genetic_params = len(relevant_muts)


    def make_seqgraph(self, colin_thres = 5):
        '''
        code amino acid differences between sequences into a matrix
        the matrix has dimensions #measurements x #observed mutations

        Parameters
        ----------
        colin_thres : int, optional
            Description
        '''
        seq_graph = []
        titer_dist = []
        weights = []

        n_params = self.genetic_params + len(self.sera) + len(self.test_strains)
        # loop over all measurements and encode the HI model as [0,1,0,1,0,0..] vector:
        # 1-> mutation present, 0 not present, same for serum and virus effects
        for (test, ref), val in self.train_titers.items():
            if not np.isnan(val) and not np.isinf(val):
                try:
                    muts = self.get_mutations(ref[0], test)
                    if muts is None:
                        continue
                    tmp = np.zeros(n_params, dtype=int) # zero vector, ones will be filled in
                    # determine branch indices on path
                    mutation_indices = np.unique([self.relevant_muts.index(mut) for mut in muts
                                                  if mut in self.relevant_muts])
                    if len(mutation_indices): tmp[mutation_indices] = 1
                    # add serum effect for heterologous viruses
                    if test!=ref[0]:
                        tmp[self.genetic_params+self.sera.index(ref)] = 1
                    # add virus effect
                    tmp[self.genetic_params+len(self.sera)+self.test_strains.index(test)] = 1
                    # append model and fit value to lists seq_graph and titer_dist
                    seq_graph.append(tmp)
                    titer_dist.append(val)
                    # for each measurment (row in the big matrix), attach weight that accounts for representation of serum
                    weights.append(1.0/(1.0 + self.serum_Kc*self.titers.measurements_per_serum[ref]))
                except:
                    import pdb; pdb.set_trace()
                    print(test, ref, "ERROR")

        # convert to numpy arrays and save product of tree graph with its transpose for future use
        self.weights = np.sqrt(weights)
        self.titer_dist =  np.array(titer_dist)*self.weights
        self.design_matrix = (np.array(seq_graph).T*self.weights).T
        if colin_thres is not None and self.genetic_params > 0:
            self.collapse_colinear_mutations(colin_thres)
        self.TgT = np.dot(self.design_matrix.T, self.design_matrix)
        print("Found", self.design_matrix.shape, "measurements x parameters")


    def collapse_colinear_mutations(self, colin_thres):
        '''
        find colinear columns of the design matrix, collapse them into clusters

        Parameters
        ----------
        colin_thres : TYPE
            Description
        '''
        TT = self.design_matrix[:,:self.genetic_params].T
        mutation_clusters = []
        n_measurements = self.design_matrix.shape[0]
        # a greedy algorithm: if column is similar to existing cluster -> merge with cluster, else -> new cluster
        for col, mut in zip(TT, self.relevant_muts):
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
        print("dimensions of old design matrix",self.design_matrix.shape)
        self.design_matrix = np.hstack((np.array([c[0] for c in mutation_clusters]).T,
                                     self.design_matrix[:,self.genetic_params:]))
        self.genetic_params = len(mutation_clusters)
        # use the first mutation of a cluster to index the effect
        # make a dictionary that maps this effect to the cluster
        self.mutation_clusters = {c[1][0]:c[1] for c in mutation_clusters}
        self.relevant_muts = [c[1][0] for c in mutation_clusters]
        print("dimensions of new design matrix",self.design_matrix.shape)


    def train(self,**kwargs):
        '''
        determine the model parameters. the result will be stored in self.substitution_effect

        Parameters
        ----------
        **kwargs
            Description
        '''
        self._train(**kwargs)
        for mi, mut in enumerate(self.relevant_muts):
            self.substitution_effect[mut] = self.model_params[mi]


    def predict_titer(self, virus, serum, cutoff=0.0):
        muts= self.get_mutations(serum[0], virus)
        if muts is not None:
            pot = self.serum_potency[serum] if serum in self.serum_potency else 0.0
            avi = self.virus_effect[virus] if virus in self.virus_effect else 0.0
            return avi + pot\
                + np.sum([self.substitution_effect[mut] for mut in muts
                if (mut in self.substitution_effect and self.substitution_effect[mut]>cutoff)])
        else:
            return None

    def compile_substitution_effects(self, cutoff=1e-4):
        '''
        compile a flat json of substitution effects for visualization, prune mutation without effect

        Parameters
        ----------
        cutoff : float, optional
            Description

        Returns
        -------
        TYPE
            Description
        '''
        return {mut[0]+':'+mut[1]:np.round(val,int(-np.log10(cutoff)))
                for mut, val in self.substitution_effect.items() if val>cutoff}

    def annotate_tree(self, tree):
        """Annotates antigenic advance attributes to nodes of a given tree built from
        the same sequences used to train the model.

        Parameters
        ----------
        tree : Bio.Phylo

        Returns
        -------
        Bio.Phylo
            input tree instance with nodes annotated by per-branch and
            cumulative antigenic advance attributes `dTiterSub` and
            `cTiterSub`
        """
        tree.root.dTiterSub = 0
        tree.root.cTiterSub = 0
        for node in tree.find_clades():
            for child in node.clades:
                # Get mutations between the current node and its parent.
                mutations = self.get_mutations(child.name, node.name)

                # Calculate titer drop on the branch to the current node.
                child.dTiterSub = 0
                for gene, mutation in mutations:
                    child.dTiterSub += self.substitution_effect.get((gene, mutation), 0)

                # Calculate the cumulative titer drop from the root to the current node.
                child.cTiterSub = node.cTiterSub + child.dTiterSub

        return tree


if __name__=="__main__":
    # test tree model (assumes there is a tree called flu in memory...)
    ttm = TreeModel(flu.tree.tree, flu.titers)
    ttm.prepare(training_fraction=0.8)
    ttm.train(method='nnl1reg')
    ttm.validate(plot=True)

    tsm = SubstitutionModel(flu.tree.tree, flu.titers)
    tsm.prepare(training_fraction=0.8)
    tsm.train(method='nnl1reg')
    tsm.validate(plot=True)
