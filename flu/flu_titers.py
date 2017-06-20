import numpy as np
from pdb import set_trace
import os

def HI_model(self):
    '''
    estimate a tree and substitution model using titers titer_fname.
    '''
    from base.titer_model import tree_model, substitution_model
    ## define the filenames here
    fname = str(self.config["titers"]["fname"].replace("<LINEAGE>", self.info["lineage"]))

    ## define the kwargs explicitly
    kwargs = {}
    if "criterium" in self.config["titers"]:
        kwargs["criterium"] = self.config["titers"]["criterium"]

    ## TREE MODEL
    self.HI_tree = tree_model(self.tree.tree, titer_fname = fname, **kwargs)
    self.HI_tree.prepare(**kwargs)
    self.HI_tree.train(**kwargs)
    # add tree attributes to the list of attributes that are saved in intermediate files
    for n in self.tree.tree.find_clades():
        n.attr['cTiter'] = n.cTiter
        n.attr['dTiter'] = n.dTiter
        # print("cumulative: {} delta: {}".format(n.cTiter, n.dTiter))
    self.config["auspice"]["color_options"]["cTiter"] = {
        "menuItem": "antigenic advance", "type": "continuous", "legendTitle": "Antigenic Advance", "key": "cTiter"
    }

    # SUBSTITUTION MODEL
    ## currently broken TODO
    # self.HI_subs = substitution_model(self.tree.tree, titer_fname = fname,**kwargs)
    # self.HI_subs.prepare(**kwargs)
    # self.HI_subs.train(**kwargs)


def HI_export(self):
    from base.io_util import write_json
    prefix = os.path.join(self.config["output"]["auspice"], self.info["prefix"])
    if hasattr(self, 'HI_tree'):
        # export the raw titers
        hi_data = self.HI_tree.compile_titers()
        write_json(hi_data, prefix+'_titers.json')
        # export the tree model (avidities and potencies only)
        tree_model = {'potency':self.HI_tree.compile_potencies(),
                      'avidity':self.HI_tree.compile_virus_effects(),
                      'dTiter':{n.clade:n.dTiter for n in self.tree.tree.find_clades() if n.dTiter>1e-6}}
        write_json(tree_model, prefix+'_titer_tree_model.json')
    else:
        print('Tree model not yet trained')

    if hasattr(self, 'HI_subs'):
        # export the substitution model
        subs_model = {'potency':self.HI_subs.compile_potencies(),
                      'avidity':self.HI_subs.compile_virus_effects(),
                      'substitution':self.HI_subs.compile_substitution_effects()}
        write_json(subs_model, prefix+'_titer_subs_model.json')
    else:
        print('Substitution model not yet trained')

def H3N2_scores(tree, epitope_mask_file, epitope_mask_version='wolf'):
    '''
    takes a H3N2 HA tree and assigns H3 specific characteristics to
    internal and external nodes
    '''
    def epitope_sites(aa):
        return aa[epitope_mask[:len(aa)]]

    def nonepitope_sites(aa):
        return aa[~epitope_mask[:len(aa)]]

    def receptor_binding_sites(aa):
        '''
        Receptor binding site mutations from Koel et al. 2014
        These are (145, 155, 156, 158, 159, 189, 193) in canonical HA numbering
        need to subtract one since python arrays start at 0
        '''
        sp = 16
        rbs = map(lambda x:x+sp-1, [145, 155, 156, 158, 159, 189, 193])
        return np.array([aa[pos] for pos in rbs])

    def get_total_peptide(node):
        '''
        the concatenation of signal peptide, HA1, HA1
        '''
        return np.fromstring(node.translations['SigPep']+node.translations['HA1']
                       + node.translations['HA2'], 'S1')

    def epitope_distance(aaA, aaB):
        """Return distance of sequences aaA and aaB by comparing epitope sites"""
        epA = epitope_sites(aaA)
        epB = epitope_sites(aaB)
        distance = np.sum(epA!=epB)
        return distance

    def nonepitope_distance(aaA, aaB):
        """Return distance of sequences aaA and aaB by comparing non-epitope sites"""
        neA = nonepitope_sites(aaA)
        neB = nonepitope_sites(aaB)
        distance = np.sum(neA!=neB)
        return distance

    def receptor_binding_distance(aaA, aaB):
        """Return distance of sequences aaA and aaB by comparing receptor binding sites"""
        neA = receptor_binding_sites(aaA)
        neB = receptor_binding_sites(aaB)
        distance = np.sum(neA!=neB)
        return distance

    epitope_map = {}
    with open(epitope_mask_file) as f:
        for line in f:
            (key, value) = line.strip().split()
            epitope_map[key] = value
    if epitope_mask_version in epitope_map:
        epitope_mask = np.fromstring(epitope_map[epitope_mask_version], 'S1')=='1'
    root = tree.root
    root_total_aa_seq = get_total_peptide(root)
    for node in tree.find_clades():
        total_aa_seq = get_total_peptide(node)
        node.attr['ep'] = epitope_distance(total_aa_seq, root_total_aa_seq)
        node.attr['ne'] = nonepitope_distance(total_aa_seq, root_total_aa_seq)
        node.attr['rb'] = receptor_binding_distance(total_aa_seq, root_total_aa_seq)
