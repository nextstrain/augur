import Bio
import time
import numpy as np
import pandas as pd
from scipy.stats import linregress
import sys

try:
    import itertools.izip as zip
except ImportError:
    pass

from .scores import calculate_LBI, select_nodes_in_season
from .titer_model import SubstitutionModel, TiterCollection, TreeModel

# all fitness predictors should be designed to give a positive sign, ie.
# number of epitope mutations
# -1 * number of non-epitope mutations

def inverse_cross_immunity_amplitude(d_ep, d_init):
    """Return the inverse cross-immunity amplitude corresponding to the given
    epitope distance between two amino acid sequences and a predetermined
    scaling parameter that controls the time period across which cross-immunity
    decays.

    Note that this implementation differs from Luksza and Lassig in that
    decaying cross-immunity is measured on a scale of 0 - 1 where no epitope
    differences correspond to a cross-immunity of 1 and more mutations decrease
    the cross-immunity score. These values are subtracted from 1 such that the
    fitness predictor in the model has positive, increasing values as
    cross-immunity wanes.
    """
    return 1 - np.exp(-d_ep / float(d_init))


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
        return "".join([value for key, value in node.translations.items() if key != "nuc"])

    def setup_predictor(self, tree, pred, timepoint, **kwargs):
        if pred == 'lbi':
            select_nodes_in_season(tree, timepoint, time_window=kwargs["time_window"])
            calculate_LBI(tree, **kwargs)

        if pred == 'cLbi':
            select_nodes_in_season(tree, timepoint, time_window=kwargs["time_window"])
            calculate_LBI(tree, attr=pred, **kwargs)

            max_lbi = 0.0
            for node in tree.find_clades():
                if node.up:
                    node.attr[pred] = node.up.attr[pred] + node.attr[pred]

                if node.attr[pred] > max_lbi:
                    max_lbi = node.attr[pred]

            for node in tree.find_clades():
                if max_lbi > 0:
                    node.attr[pred] = node.attr[pred] / max_lbi

                setattr(node, pred, node.attr[pred])

        if pred == 'ep':
            self.calc_epitope_distance(tree)
        if pred == 'ep_x':
            self.calc_epitope_cross_immunity(tree, timepoint, **kwargs)
        #if pred == 'ne':
        #    self.calc_nonepitope_distance(tree)
        if pred == 'ne_star':
            self.calc_nonepitope_star_distance(tree, timepoint, **kwargs)
        if pred == 'tol':
            self.calc_tolerance(tree, preferences_file='metadata/2017-12-07-H3N2-preferences-rescaled.csv', attr = 'tol')
        if pred == 'tol_mask':
            self.calc_tolerance(tree, preferences_file='metadata/2017-12-07-H3N2-preferences-rescaled.csv', attr = pred, use_epitope_mask=True)
        if pred == 'dms':
            if not hasattr(tree.root, pred):
                self.calc_dms(tree, kwargs.get("preferences_file"))
        if pred == 'tol_ne':
            self.calc_tolerance(tree, epitope_mask = self.tolerance_mask, attr = 'tol_ne')
        if pred == 'null':
            self.calc_null_predictor(tree)
        if pred == 'random':
            self.calc_random_predictor(tree)
        #if pred == 'dfreq':
            # do nothing
        if pred == 'cTiter':
            self.calc_titer_model("tree", tree, timepoint, **kwargs)
        if pred == 'cTiterSub':
            self.calc_titer_model("substitution", tree, timepoint, **kwargs)
        if pred == 'future_fitness':
            self.calc_future_fitness(tree, timepoint, **kwargs)

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
        assert len(aa) == len(self.epitope_mask), "Sequence length: %s long, mask length: %s" % (len(aa), len(self.epitope_mask))
        sites = []
        for a, m in zip(aa, self.epitope_mask):
            if m == '1':
                sites.append(a)
        return ''.join(sites)

    def nonepitope_sites(self, aa):
        """Returns amino acids from the given protein sequence corresponding to
        non-epitope sites.
        """
        sites = []
        for a, m in zip(aa, self.tolerance_mask):
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
        distance = sum(a != b for a, b in zip(epA, epB))
        return distance

    def fast_epitope_distance(self, epA, epB):
        """Return distance of sequences aaA and aaB by comparing epitope sites"""
        return np.count_nonzero(epA!=epB)

    def nonepitope_distance(self, aaA, aaB):
        """Return distance of sequences aaA and aaB by comparing non-epitope sites"""
        neA = self.nonepitope_sites(aaA)
        neB = self.nonepitope_sites(aaB)
        distance = sum(a != b for a, b in zip(neA, neB))
        return distance

    def rbs_distance(self, aaA, aaB):
        """Return distance of sequences aaA and aaB by comparing receptor binding sites (Koel sites)"""
        rbsA = self.receptor_binding_sites(aaA)
        rbsB = self.receptor_binding_sites(aaB)
        distance = sum(a != b for a, b in zip(rbsA, rbsB))
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

    def calc_epitope_cross_immunity(self, tree, timepoint, step_size, d_init=14, attr='ep_x', years_to_wane=5, **kwargs):
        """Calculates the distance at epitope sites to contemporaneous viruses
        this should capture cross-immunity of circulating viruses
        meant to be used in conjunction with epitope_distance that focuses
        on escape from previous human immunity.
        """
        # Find all strains sampled between the previous timepoint and the
        # current timepoint. We will calculate cross-immunity for each of these
        # strains. Additionally, find all strains sampled prior to the previous
        # timepoint. We will compare the current timepoint's strains to these
        # past strains.
        previous_timepoint = timepoint - step_size
        current_nodes = []
        past_nodes = []
        for node in tree.get_terminals():
            if node.attr['num_date'] < previous_timepoint:
                past_nodes.append(node)
            elif previous_timepoint <= node.attr['num_date'] < timepoint:
                current_nodes.append(node)

            # Annotate an array of amino acids at epitope sites to each node.
            if not hasattr(node, 'np_ep'):
                if not hasattr(node, 'aa'):
                    node.aa = self._translate(node)
                node.np_ep = np.array(list(self.epitope_sites(node.aa)))

            # Initialize all tips to 0 distance.
            setattr(node, attr, 0.0)

        print("calculating cross-immunity at %s between %s current nodes and %s past nodes" % (timepoint, len(current_nodes), len(past_nodes)))

        for node in current_nodes:
            waning_effects = []
            frequencies = []
            distances = []
            for comp_node in past_nodes:
                # Consider only past viruses that were sampled while immunity
                # had not waned relative to the current virus.
                if node.attr["num_date"] - comp_node.attr["num_date"] < years_to_wane:
                    # Calculate the effect of linear waning immunity on the
                    # overall cross-immunity amplitude as a proportion of the
                    # years between the two viruses and the overall waning
                    # period.
                    waning_effect = 1 - ((node.attr["num_date"] - comp_node.attr["num_date"]) / years_to_wane)
                    waning_effects.append(waning_effect)

                    # Track the frequency of each past node and its distance from the current node.
                    # Cross-immunity is scaled by the maximum frequency that the past node ever obtained.
                    frequencies.append(max(comp_node.censored_freqs.values()))

                    # Calculate the number of epitope mutations between viruses.
                    distances.append(self.fast_epitope_distance(node.np_ep, comp_node.np_ep))

            # Calculate inverse cross-immunity amplitude once from all distances to the current strain.
            # This is an increasingly positive value for strains that are increasingly distant from previous strains.
            cross_immunity_amplitudes = inverse_cross_immunity_amplitude(np.array(distances), d_init)

            # Scale cross-immunity by waning effects and past strain frequencies
            # and sum across all past viruses.
            waning_effects = np.array(waning_effects)
            frequencies = np.array(frequencies)
            total_cross_immunity = (waning_effects * frequencies * cross_immunity_amplitudes).sum()
            setattr(node, attr, total_cross_immunity)

    def calc_nonepitope_star_distance(self, tree, timepoint, step_size, attr='ne_star', **kwargs):
        """Calculate the non-epitope mutation distance between each node in the current
        timepoint interval and its closest ancestor in the previous timepoint
        interval.
        """
        # First, find all nodes sampled in this timepoint interval and the
        # previous one and annotate non-epitope mutations to each node.
        previous_timepoint = timepoint - step_size
        current_nodes = []
        for node in tree.find_clades():
            if not hasattr(node, 'np_ne'):
                if not hasattr(node, 'aa'):
                    node.aa = self._translate(node)
                node.np_ne = np.array(list(self.nonepitope_sites(node.aa)))

            if node.is_terminal():
                # Initialize predictor to zero.
                if not hasattr(node, attr):
                    setattr(node, attr, 0.0)

                if previous_timepoint <= node.attr['num_date'] < timepoint:
                    current_nodes.append(node)

        # Next, find the first ancestor of each current node that was sampled in
        # the previous time interval.
        for node in current_nodes:
            parent = node.up
            while parent.attr["num_date"] > previous_timepoint:
                parent = parent.up

            # Calculate the non-epitope distance to the ancestor.
            setattr(node, attr, self.fast_epitope_distance(node.np_ne, parent.np_ne))

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

    def calc_titer_model(self, model, tree, timepoint, titers, lam_avi, lam_pot, lam_drop, **kwargs):
        """Calculates the requested titer model for the given tree using only titers
        associated with strains sampled prior to the given timepoint.
        """
        # Filter titers by date from the root node to the current timepoint.
        node_lookup = {node.name: node for node in tree.get_terminals()}
        filtered_titers = TiterCollection.subset_to_date(titers, node_lookup, tree.root.numdate, timepoint)

        # Set titer model parameters.
        kwargs = {
            "criterium": lambda node: hasattr(node, "aa_muts") and sum([len(node.aa_muts[protein]) for protein in node.aa_muts]) > 0,
            "lam_avi": lam_avi,
            "lam_pot": lam_pot,
            "lam_drop": lam_drop
        }

        # Run the requested model and annotate results to nodes.
        if model == "tree":
            model = TreeModel(tree, filtered_titers, **kwargs)
        elif model == "substitution":
            model = SubstitutionModel(tree, filtered_titers, **kwargs)

        model.prepare(**kwargs)
        model.train(**kwargs)

    def calc_future_fitness(self, tree, timepoint, attr="future_fitness", **kwargs):
        """Calculate the known future frequency of each tip at the given timepoint.

        This predictor is a positive control for the model that should always
        predict the correct future frequencies since it is borrowing that
        information from the future without any censoring.
        """
        for node in tree.get_terminals():
            # Try to use the known future frequency of a node and fallback to
            # the current timepoint frequency when we don't know the
            # future. This should only be true for the last timepoint.
            if timepoint in node.observed_final_freqs:
                setattr(node, attr, node.observed_final_freqs[timepoint])
            else:
                setattr(node, attr, node.timepoint_freqs[timepoint])
