import numpy as np
import re
import sys


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

def calculate_metadata_scores(tree):
    """Calculate scores for each node in the given tree based on metadata already
    annotated on each node including age and sex of infected individuals.
    """
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

def mask_sites(aa, mask):
    return aa[mask[:len(aa)]]

def mask_distance(aaA, aaB, mask):
    """Return distance of sequences aaA and aaB by comparing sites in the given binary mask."""
    sites_A = mask_sites(aaA, mask)
    sites_B = mask_sites(aaB, mask)
    distance = np.sum(sites_A != sites_B)
    return distance

def glycosylation_count(total_aa_seq, glyc_mask):
    # TODO: need to restrict to surface residues.
    total_aa_seq_masked = "".join([aa if mask else 'X'
                                   for (mask, aa) in zip(glyc_mask, total_aa_seq)])

    return len(re.findall('N[^P][ST][^P]', total_aa_seq_masked))

def calculate_sequence_scores(tree, mask_file, lineage, segment, epitope_mask_version='wolf', glyc_mask_version='wolf'):
    """Calculate scores from the amino acid sequence of each node in the given tree.

    Sequence scores depend on lineage- and segment-specific amino acid site masks or named masks.
    """
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

    # Annotate scores to each node based on its amino acid sequence.
    for node in tree.find_clades():
        total_aa_seq = get_total_peptide(node, segment)
        node.attr['ep'] = mask_distance(total_aa_seq, root_total_aa_seq, epitope_mask)
        node.attr['ne'] = mask_distance(total_aa_seq, root_total_aa_seq, ~epitope_mask)
        node.attr['glyc'] = glycosylation_count(total_aa_seq, glyc_mask)

        if rbs_mask is not None:
            node.attr['rb'] = mask_distance(total_aa_seq, root_total_aa_seq, rbs_mask)
