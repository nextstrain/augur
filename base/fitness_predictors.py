import Bio
import time
import numpy as np
import pandas as pd
from itertools import izip
from scipy.stats import linregress
import sys

from builds.flu.scores import calculate_LBI

# all fitness predictors should be designed to give a positive sign, ie.
# number of epitope mutations
# -1 * number of non-epitope mutations

class fitness_predictors(object):

    def __init__(self, predictor_names = ['ep', 'lb', 'dfreq'], **kwargs):
        if "epitope_masks_fname" in kwargs and "epitope_mask_version" in kwargs and "tolerance_mask_version" in kwargs:
            self.setup_epitope_mask(epitope_masks_fname = kwargs["epitope_masks_fname"], epitope_mask_version = kwargs["epitope_mask_version"], tolerance_mask_version = kwargs["tolerance_mask_version"])
        elif "epitope_masks_fname" in kwargs and "epitope_mask_version" in kwargs:
            self.setup_epitope_mask(epitope_masks_fname = kwargs["epitope_masks_fname"], epitope_mask_version = kwargs["epitope_mask_version"])
        else:
            self.setup_epitope_mask()
        self.predictor_names = predictor_names

    def _translate(self, node):
        """Returns a single amino acid sequence corresponding to the given node's
        protein translations concatenated in genomic order (e.g., SigPep, HA1,
        and HA2 for flu's HA protein).
        """
        return "".join(node.translations.values())

    def setup_predictor(self, tree, pred, timepoint, **kwargs):
        if pred == 'lb':
           calculate_LBI(tree, **kwargs)
        # if pred == 'ep':
        #     self.calc_epitope_distance(tree)
        if pred == 'ep_x':
            self.calc_epitope_cross_immunity(tree, timepoint)
        #if pred == 'ne':
        #    self.calc_nonepitope_distance(tree)
        if pred == 'ne_star':
            self.calc_nonepitope_star_distance(tree)
        if pred == 'tol':
            self.calc_tolerance(tree, preferences_file='metadata/2017-12-07-H3N2-preferences-rescaled.csv', attr = 'tol')
        if pred == 'tol_mask':
            self.calc_tolerance(tree, preferences_file='metadata/2017-12-07-H3N2-preferences-rescaled.csv', attr = pred, use_epitope_mask=True)
        if pred == 'dms':
            if not hasattr(tree.root, pred):
                self.calc_dms(tree, preferences_file='metadata/2017-12-07-H3N2-preferences-rescaled.csv')
        if pred == 'tol_ne':
            self.calc_tolerance(tree, epitope_mask = self.tolerance_mask, attr = 'tol_ne')
        if pred == 'null':
            self.calc_null_predictor(tree)
        if pred == 'random':
            self.calc_random_predictor(tree)
        #if pred == 'dfreq':
            # do nothing
        #if pred == 'cHI':
            # do nothing

    def setup_epitope_mask(self, epitope_masks_fname = 'metadata/ha_masks.tsv', epitope_mask_version = 'wolf', tolerance_mask_version = 'ha1'):
        sys.stderr.write("setup " + str(epitope_mask_version) + " epitope mask and " + str(tolerance_mask_version) + " tolerance mask\n")
        self.epitope_mask = ""
        self.tolerance_mask = ""
        epitope_map = {}
        with open(epitope_masks_fname) as f:
            for line in f:
                (key, value) = line.split()
                epitope_map[key] = value
        if epitope_mask_version in epitope_map:
            self.epitope_mask = epitope_map[epitope_mask_version]
        if tolerance_mask_version in epitope_map:
            self.tolerance_mask = epitope_map[tolerance_mask_version]
        else:
            self.tolerance_mask = self.epitope_mask

    def epitope_sites(self, aa):
        """Returns amino acids from the given protein sequence corresponding to sites in
        a predefined epitope mask.
        """
        sites = []
        for a, m in izip(aa, self.epitope_mask):
            if m == '1':
                sites.append(a)
        return ''.join(sites)

    def nonepitope_sites(self, aa):
        """Returns amino acids from the given protein sequence corresponding to
        non-epitope sites.
        """
        sites = []
        for a, m in izip(aa, self.tolerance_mask):
            if m == '0':
                sites.append(a)
        return ''.join(sites)

    def receptor_binding_sites(self, aa):
        """Returns amino acids from the given protein sequence corresponding to seven
        Koel et al. 2013 receptor binding sites.
        """
        sp = 16
        aaa = np.fromstring(aa, 'S1')
        receptor_binding_list = map(lambda x:x+sp-1, [145, 155, 156, 158, 159, 189, 193])
        return ''.join(aaa[receptor_binding_list])

    def epitope_distance(self, aaA, aaB):
        """Return distance of sequences aaA and aaB by comparing epitope sites"""
        epA = self.epitope_sites(aaA)
        epB = self.epitope_sites(aaB)
        distance = sum(a != b for a, b in izip(epA, epB))
        return distance

    def fast_epitope_distance(self, epA, epB):
        """Return distance of sequences aaA and aaB by comparing epitope sites"""
        return np.count_nonzero(epA!=epB)

    def nonepitope_distance(self, aaA, aaB):
        """Return distance of sequences aaA and aaB by comparing non-epitope sites"""
        neA = self.nonepitope_sites(aaA)
        neB = self.nonepitope_sites(aaB)
        distance = sum(a != b for a, b in izip(neA, neB))
        return distance

    def rbs_distance(self, aaA, aaB):
        """Return distance of sequences aaA and aaB by comparing receptor binding sites (Koel sites)"""
        rbsA = self.receptor_binding_sites(aaA)
        rbsB = self.receptor_binding_sites(aaB)
        distance = sum(a != b for a, b in izip(rbsA, rbsB))
        return distance

    def calc_epitope_distance(self, tree, attr='ep', ref = None):
        '''
        calculates the distance at epitope sites of any tree node to ref
        tree   --   dendropy tree
        attr   --   the attribute name used to save the result
        '''
        for node in tree.find_clades(order="postorder"):
            if not hasattr(node, 'np_ep'):
                if not hasattr(node, 'aa'):
                    node.aa = self._translate(node)
                node.np_ep = np.array(list(self.epitope_sites(node.aa)))
        if ref == None:
            ref = tree.root
        for node in tree.find_clades(order="postorder"):
            node.__setattr__(attr, self.fast_epitope_distance(node.np_ep, ref.np_ep))

    def calc_epitope_cross_immunity(self, tree, timepoint, window = 2.0, attr='ep_x'):
        '''
        calculates the distance at epitope sites to contemporaneous viruses
        this should capture cross-immunity of circulating viruses
        meant to be used in conjunction with epitope_distance that focuses
        on escape from previous human immunity
        tree   --   dendropy tree
        attr   --   the attribute name used to save the result
        '''
        comparison_nodes = []
        for node in tree.find_clades(order="postorder"):
            if node.is_leaf():
                if node.attr['num_date'] < timepoint and node.attr['num_date'] > timepoint - window:
                    comparison_nodes.append(node)
            if not hasattr(node, 'np_ep'):
                if not hasattr(node, 'aa'):
                    node.aa = self._translate(node)
                node.np_ep = np.array(list(self.epitope_sites(node.aa)))
        print "calculating cross-immunity to " + str(len(comparison_nodes)) + " comparison nodes"
        for node in tree.find_clades(order="postorder"):
            mean_distance = 0
            count = 0
            for comp_node in comparison_nodes:
                mean_distance += self.fast_epitope_distance(node.np_ep, comp_node.np_ep)
                count += 1
            if count > 0:
                mean_distance /= float(count)
            node.__setattr__(attr, mean_distance)

    def calc_rbs_distance(self, tree, attr='rb', ref = None):
        '''
        calculates the distance at receptor binding sites of any tree node to ref
        tree   --   dendropy tree
        attr   --   the attribute name used to save the result
        '''
        if ref == None:
            ref = self._translate(tree.root)
        for node in tree.find_clades(order="postorder"):
            if not hasattr(node, 'aa'):
                node.aa = self._translate(node)
            node.__setattr__(attr, self.rbs_distance(node.aa, ref))

    def calc_dms(self, tree, preferences_file, attr="dms"):
        """Calculates the cumulative mutational effect of mutations per node relative to
        their parent node in the given tree and based on a set of site-specific
        AA preferences.

        The calculation can be limited to a given list of `positions` in the
        given sequence.
        """
        preferences = pd.read_csv(preferences_file)
        preferences = preferences.reset_index()
        stacked_preferences = preferences.loc[:, "A":"Y"].stack()

        # Use all positions in the given sequence by default.
        tree.root.aa = self._translate(tree.root)
        positions = list(range(len(tree.root.aa)))

        # Calculate a default value for missing preferences.
        missing_preference = 1e-10

        for node in tree.root.find_clades():
            # Determine amino acid sequence if it is not defined.
            if not hasattr(node, "aa"):
                node.aa = self._translate(node)

            # Get amino acid mutations between this node and its parent.
            if node.up is None:
                mut_effect = 0.0
            else:
                parent_mutations = []
                node_mutations = []

                for i in positions:
                    if node.aa[i] != node.up.aa[i]:
                        parent_mutations.append((i, node.up.aa[i]))
                        node_mutations.append((i, node.aa[i]))

                # Calculate mutational effect for this node based on its differences
                # since the parent node. If there are no mutations since the parent,
                # this node keeps the same effect. Otherwise, sum the effects of all
                # mutations since the parent.
                mut_effect = node.up.attr[attr]
                if len(node_mutations) > 0:
                    # Mutational effect is the preference of the current node's amino acid
                    # at each mutated site divided by the preference of the original and
                    # transformed to log scale such that changes to less preferred amino acids
                    # will produce a negative effect and changes to more preferred will be positive.
                    node_preferences = stacked_preferences[node_mutations].fillna(missing_preference)
                    parent_preferences = stacked_preferences[parent_mutations].fillna(missing_preference)
                    new_mut_effect = np.log2(node_preferences.values / parent_preferences.values).sum()
                    mut_effect += new_mut_effect

            setattr(node, attr, mut_effect)
            node.attr[attr] = mut_effect

    def calc_tolerance(self, tree, preferences_file, attr="tol", use_epitope_mask=False):
        """Calculates log odds of a node's AA sequence relative to a set of
        site-specific AA preferences.

        Takes a tree and a path to a DMS preferences file with one row per site
        and labeled columns per amino acid.
        """
        preferences = pd.read_csv(preferences_file)
        preferences = preferences.reset_index()
        stacked_preferences = preferences.loc[:, "A":"Y"].stack()

        # Use all positions in the given sequence by default.
        tree.root.aa = self._translate(tree.root)
        if use_epitope_mask:
            positions = [i for i in range(len(self.epitope_mask))
                         if self.epitope_mask[i] == "1"]
            print("Using %i positions for sequence preferences from epitope mask." % (len(positions)))
        else:
            positions = list(range(len(tree.root.aa)))
            print("Using %i positions for sequence preferences from all sites." % (len(positions)))

        # Log scale all preferences prior to tolerance calculation.
        log_preferences = np.log(stacked_preferences)

        # Calculate a default value for missing preferences.
        missing_preference = np.log(1e-10)

        # Calculate the sequence preference for the root node.

        # Create a list of keys into the preference series where the first
        # value is the zero-based protein position and the second value is
        # the amino acid at that position in the given sequence.
        aa_array = [(i, tree.root.aa[i]) for i in positions]

        # Look up preferences for the amino acids at each site in the given
        # protein using the given `preferences` series indexed by site and
        # amino acid.
        node_preferences = log_preferences[aa_array]

        # Replace missing values with a very small probability. This
        # primarily accounts for stop codons ("X") that do not have an
        # associated DMS preference.
        node_preferences = node_preferences.fillna(missing_preference)

        # Calculate sum of the log of the preferences the node's amino acid
        # sequence.
        tolerance = node_preferences.sum() / len(positions)
        setattr(tree.root, attr, tolerance)
        tree.root.attr[attr] = tolerance

        # Calculate the tolerance for the remaining nodes in the tree based
        # on the difference between each node and its parent.
        for node in tree.root.find_clades():
            # Determine amino acid sequence if it is not defined.
            if not hasattr(node, "aa"):
                node.aa = self._translate(node)

            # Skip the root.
            if node == tree.root:
                continue

            # Get amino acid mutations between this node and its parent.
            node_mutations = []
            parent_mutations = []
            for i in positions:
                if node.aa[i] != node.up.aa[i]:
                    parent_mutations.append((i, node.up.aa[i]))
                    node_mutations.append((i, node.aa[i]))

            # Assign a new tolerance to this node based on its differences
            # since the parent node.
            tolerance = node.up.attr[attr]
            if len(node_mutations) > 0:
                tolerance = (tolerance -
                             (log_preferences[parent_mutations].fillna(missing_preference).sum() / len(positions)) +
                             (log_preferences[node_mutations].fillna(missing_preference).sum() / len(positions)))

            setattr(node, attr, tolerance)
            node.attr[attr] = tolerance

    def calc_nonepitope_distance(self, tree, attr='ne', ref = None):
        '''
        calculates -1 * distance at nonepitope sites of each node in tree to ref
        tree   --   dendropy tree
        attr   --   the attribute name used to save the result
        '''
        if ref == None:
            ref = self._translate(tree.root)
        for node in tree.find_clades(order="postorder"):
            if not hasattr(node, 'aa'):
                node.aa = self._translate(node)
            distance = self.nonepitope_distance(node.aa, ref)
            node.__setattr__(attr, -1 * distance)

    def calc_nonepitope_star_distance(self, tree, attr='ne_star', seasons = []):
        '''
        calculates the distance at nonepitope sites of any tree node to ref
        tree   --   dendropy tree
        attr   --   the attribute name used to save the result
        '''
        for node in tree.find_clades(order="postorder"):
            if len(node.season_tips) and node!=tree.root:
                if not hasattr(node, 'aa'):
                    node.aa = self._translate(node)
                tmp_node = node.parent_node
                cur_season = min(node.season_tips.keys())
                prev_season = seasons[max(0,seasons.index(cur_season)-1)]
                while True:
                    if tmp_node!=tree.root:
                        if prev_season in tmp_node.season_tips and len(tmp_node.season_tips[prev_season])>0:
                            break
                        else:
                            tmp_node=tmp_node.parent_node
                    else:
                        break
                if not hasattr(tmp_node, 'aa'):
                    tmp_node.aa = self._translate(tmp_node)
                node.__setattr__(attr, self.nonepitope_distance(node.aa, tmp_node.aa))
            else:
                node.__setattr__(attr, np.nan)

    def calc_null_predictor(self, tree, attr="null"):
        """Assign a zero value to each node as a control representing a null predictor.
        """
        for node in tree.find_clades():
            setattr(node, attr, 0)

    def calc_random_predictor(self, tree, attr="random"):
        """Assign a random value between 0 and 1 to each node as a control predictor.
        """
        for node in tree.find_clades():
            setattr(node, attr, np.random.random())
