import numpy as np
from pdb import set_trace
import os,re
import sys

def HI_model(process):
    '''
    estimate a tree and substitution model using titers.
    '''
    from base.titer_model import TreeModel, SubstitutionModel

    ## define the kwargs explicitly
    kwargs = process.config["titers"]
    # if "criterium" in process.config["titers"]:
    #     kwargs["criterium"] = process.config["titers"]["criterium"]

    ## TREE MODEL
    process.HI_tree = TreeModel(process.tree.tree, process.titers, **kwargs)
    process.HI_tree.prepare(**kwargs)
    process.HI_tree.train(**kwargs)
    # add tree attributes to the list of attributes that are saved in intermediate files
    for n in process.tree.tree.find_clades():
        n.attr['cTiter'] = n.cTiter
        n.attr['dTiter'] = n.dTiter
        # print("cumulative: {} delta: {}".format(n.cTiter, n.dTiter))
    process.config["auspice"]["color_options"]["cTiter"] = {
        "menuItem": "antigenic advance", "type": "continuous", "legendTitle": "Antigenic Advance", "key": "cTiter"
    }

    # SUBSTITUTION MODEL
    process.HI_subs = SubstitutionModel(process.tree.tree, process.titers, **kwargs)
    process.HI_subs.prepare(**kwargs)
    process.HI_subs.train(**kwargs)

    for node in process.tree.tree.find_clades():
        dTiterSub = 0

        if hasattr(node, "aa_muts"):
            for gene, mutations in node.aa_muts.iteritems():
                for mutation in mutations:
                    dTiterSub += process.HI_subs.substitution_effect.get((gene, mutation), 0)

        node.dTiterSub = dTiterSub

        if node.up is not None:
            node.cTiterSub = node.up.cTiterSub + dTiterSub
        else:
            node.cTiterSub = 0

        node.attr["cTiterSub"] = node.cTiterSub
        node.attr["dTiterSub"] = node.dTiterSub


def HI_export(process):
    from base.io_util import write_json
    prefix = os.path.join(process.config["output"]["auspice"], process.info["prefix"])
    if hasattr(process, 'HI_tree'):
        # export the raw titers
        hi_data = process.HI_tree.compile_titers()
        write_json(hi_data, prefix+'_titers.json')
        # export the tree model (avidities and potencies only)
        tree_model = {'potency':process.HI_tree.compile_potencies(),
                      'avidity':process.HI_tree.compile_virus_effects(),
                      'dTiter':{n.clade:n.dTiter for n in process.tree.tree.find_clades() if n.dTiter>1e-6}}
        write_json(tree_model, prefix+'_titer_tree_model.json')
    else:
        print('Tree model not yet trained')

    if hasattr(process, 'HI_subs'):
        # export the substitution model
        subs_model = {'potency':process.HI_subs.compile_potencies(),
                      'avidity':process.HI_subs.compile_virus_effects(),
                      'substitution':process.HI_subs.compile_substitution_effects()}
        write_json(subs_model, prefix+'_titer_subs_model.json')
    else:
        print('Substitution model not yet trained')


def get_total_peptide(node, segment='ha'):
    '''
    the concatenation of signal peptide, HA1, HA1
    '''
    if segment=='ha':
        return np.fromstring(node.translations['SigPep'] + node.translations['HA1'] + node.translations['HA2'], 'S1')
    elif segment=='na':
        return np.fromstring(node.translations['NA'], 'S1')

def read_masks(mask_file):
    ha_masks = {}
    with open(mask_file) as f:
        for line in f:
            (key, value) = line.strip().split()
            ha_masks[key] = np.fromstring(value, 'S1')=='1'
    return ha_masks

def calculate_metadata_scores(tree, segment):
    root = tree.root

    for node in tree.get_terminals():
        node.tip_count=1.0
        if ('age' in node.attr) and type(node.attr['age']) in [str, unicode]:
            tmp_age = node.attr['age']
            if tmp_age[-1] in ['y', 'm']:
                node.attr['age'] = int(tmp_age[:-1])
                if tmp_age[-1]=='m':
                    node.attr['age']/=12.0
                elif node.attr['age']>120:
                    node.attr['age'] = 'unknown'
            else:
                node.attr['age'] = 'unknown'
        else:
            node.attr['age'] = 'unknown'
        node.tip_ages = [] if node.attr['age']=="unknown" else [node.attr['age']]

    # assign a list of ages to each internal node (will be used to calculate changes in age distribution)
    for node in tree.get_nonterminals(order='postorder'):
        valid_ages = [c.attr['age'] for c in node if c.attr['age'] is not "unknown"]
        node.tip_ages = sum( [c.tip_ages for c in node] , [])
        valid_tipcounts = [c.tip_count for c in node if c.attr['age'] is not "unknown"]
        node.tip_count = np.sum(valid_tipcounts)
        if len(valid_ages):
            node.attr['age'] = np.sum([a*t for a,t in zip(valid_ages, valid_tipcounts)])/np.sum(valid_tipcounts)
        else:
            node.attr['age'] = "unknown"
    root.outgroup_ages = []

    # assign avg_age. This is redundant in that deep nodes have average ages already
    # but terminal nodes or small clades don't have a smooth average. this assigns those
    # to the parent.
    # outgroup_ages: these are ages of all nodes NOT within the clade.
    root.attr['avg_age'] = root.attr['age']
    for node in tree.get_nonterminals(order='preorder'):
        for c1 in node:
            c1.outgroup_ages = node.outgroup_ages + sum([c2.tip_ages for c2 in node if c1!=c2], [])
            if c1.tip_count<20:
                c1.attr['avg_age'] = node.attr['avg_age']
            else:
                c1.attr['avg_age'] = c1.attr['age']

    # calculate the a score to spot changing age distributions (for what ever reason, i.e. geography)
    # these are calculated as a 2-sample KS test based on ages within and outside the clade in question
    from scipy import stats
    for node in tree.find_clades(order='preorder'):
        if node==root:
            node.attr['age_score'] = 0
            continue
        if len(node.tip_ages)>20 and len(node.outgroup_ages)>20:
            ks = stats.ks_2samp(node.outgroup_ages, node.tip_ages)
            node.attr['age_score'] = -np.log10(ks.pvalue)
        else:
            node.attr['age_score'] = node.up.attr['age_score']

        if np.isnan(node.attr['age_score']):
            node.attr['age_score']=0.0

    # gender seems not to be a relevant quantity as of now...
    for node in tree.get_terminals():
        node.tip_count=1.0
        if ('gender' in node.attr) and type(node.attr['gender']) in [str, unicode]:
            g = node.attr["gender"]
            node.attr["num_gender"] = -1 if g=='male' else (1 if g=='female' else 0)
        else:
            node.attr["num_gender"] = 0

    for node in tree.get_nonterminals(order='postorder'):
        node.tip_count = np.sum([c.tip_count for c in node])
        node.attr['num_gender'] = np.sum([c.attr['num_gender']*c.tip_count for c in node])/node.tip_count


def epitope_sites(aa, epitope_mask):
    return aa[epitope_mask[:len(aa)]]

def nonepitope_sites(aa, epitope_mask):
    return aa[~epitope_mask[:len(aa)]]

def epitope_distance(aaA, aaB, epitope_mask):
    """Return distance of sequences aaA and aaB by comparing epitope sites"""
    epA = epitope_sites(aaA, epitope_mask)
    epB = epitope_sites(aaB, epitope_mask)
    distance = np.sum(epA!=epB)
    return distance

def nonepitope_distance(aaA, aaB, epitope_mask):
    """Return distance of sequences aaA and aaB by comparing non-epitope sites"""
    neA = nonepitope_sites(aaA, epitope_mask)
    neB = nonepitope_sites(aaB, epitope_mask)
    distance = np.sum(neA!=neB)
    return distance

def glycosylation_count(total_aa_seq, glyc_mask):
    # TODO: need to restrict to surface residues.
    total_aa_seq_masked = "".join([aa if mask else 'X'
                                   for (mask, aa) in zip(glyc_mask, total_aa_seq)])

    return len(re.findall('N[^P][ST][^P]', total_aa_seq_masked))

def calculate_sequence_scores(tree, mask_file, lineage, segment, epitope_mask_version='wolf', glyc_mask_version='wolf'):
    '''
    takes a H3N2 HA tree and assigns H3 specific characteristics to
    internal and external nodes
    '''
    ha_masks = read_masks(mask_file)
    try:
        epitope_mask = ha_masks[epitope_mask_version]
    except KeyError, e:
        sys.stderr.write("ERROR: Could not find an epitope mask named '%s'.\n" % epitope_mask_version)
        raise e

    # Get amino acid sequence of root node.
    root = tree.root
    root_total_aa_seq = get_total_peptide(root, segment)

    # Setup receptor binding site mask.
    rbs_mask_name = "%s_%s_rbs" % (lineage, segment)
    rbs_mask = ha_masks.get(rbs_mask_name)

    # Setup glycosylation mask.
    if glyc_mask_version in ha_masks:
        glyc_mask = ha_masks[glyc_mask_version]
    else:
        sys.stderr.write("WARNING: Could not find a glycosylation mask named '%s'; using all positions instead.\n" % glyc_mask_version)
        glyc_mask = np.ones(len(root_total_aa_seq), dtype='bool')

    for node in tree.find_clades():
        total_aa_seq = get_total_peptide(node, segment)
        node.attr['ep'] = epitope_distance(total_aa_seq, root_total_aa_seq, epitope_mask)
        node.attr['ne'] = nonepitope_distance(total_aa_seq, root_total_aa_seq, epitope_mask)
        node.attr['glyc'] = glycosylation_count(total_aa_seq, glyc_mask)

        if rbs_mask is not None:
            node.attr['rb'] = epitope_distance(total_aa_seq, root_total_aa_seq, rbs_mask)
