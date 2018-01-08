import Bio
import time
import numpy as np
from itertools import izip
from scipy.stats import linregress
import sys

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

    def setup_predictor(self, tree, pred, timepoint):
        if pred == 'lb':
            self.calc_LBI(tree, tau = 0.0005, transform = lambda x:x)
        if pred == 'ep':
            self.calc_epitope_distance(tree)
        if pred == 'ep_x':
            self.calc_epitope_cross_immunity(tree, timepoint)
        if pred == 'ne':
            self.calc_nonepitope_distance(tree)
        if pred == 'ne_star':
            self.calc_nonepitope_star_distance(tree)
        if pred == 'tol':
            self.calc_tolerance(tree, attr = 'tol')
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

    def setup_epitope_mask(self, epitope_masks_fname = 'builds/flu/metadata/h3n2_epitope_masks.tsv', epitope_mask_version = 'wolf', tolerance_mask_version = 'ha1'):
        sys.stderr.write("setup " + epitope_mask_version + " epitope mask and " + tolerance_mask_version + " tolerance mask\n")
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
                if node.num_date < timepoint and node.num_date > timepoint - window:
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

    def calc_tolerance(self, tree, epitope_mask=None, attr='tol'):
        '''
        calculates log odds of a node's AA sequence relative to a set of site-specific AA preferences
        tree   --   dendropy tree
        attr   --   the attribute name used to save the result
        '''
        assign_fitness_tolerance(tree, epitope_mask=epitope_mask, attr=attr)

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

    def calc_LBI(self, tree, attr = 'lb', tau=0.0005, transform = lambda x:x):
        '''
        traverses the tree in postorder and preorder to calculate the
        up and downstream tree length exponentially weighted by distance.
        then adds them as LBI
        tree -- dendropy tree for whose node the LBI is being computed
        attr     -- the attribute name used to store the result
        '''
        # traverse the tree in postorder (children first) to calculate msg to parents
        for node in tree.find_clades(order="postorder"):
            node.down_polarizer = 0
            node.up_polarizer = 0
            for child in node.child_nodes():
                node.up_polarizer += child.up_polarizer
            bl =  node.edge_length/tau
            node.up_polarizer *= np.exp(-bl)
            if node.alive: node.up_polarizer += tau*(1-np.exp(-bl))

        # traverse the tree in preorder (parents first) to calculate msg to children
        for node in tree.preorder_internal_node_iter():
            for child1 in node.child_nodes():
                child1.down_polarizer = node.down_polarizer
                for child2 in node.child_nodes():
                    if child1!=child2:
                        child1.down_polarizer += child2.up_polarizer

                bl =  child1.edge_length/tau
                child1.down_polarizer *= np.exp(-bl)
                if child1.alive: child1.down_polarizer += tau*(1-np.exp(-bl))

        # go over all nodes and calculate the LBI (can be done in any order)
        for node in tree.find_clades(order="postorder"):
            tmp_LBI = node.down_polarizer
            for child in node.child_nodes():
                tmp_LBI += child.up_polarizer
            node.__setattr__(attr, transform(tmp_LBI))

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
