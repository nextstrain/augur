"""
Infer ancestral traits based on a tree.
"""

import numpy as np
from collections import defaultdict
import os, sys
import pandas as pd
from .utils import read_metadata, write_json, get_json_name
TINY = 1e-12

def mugration_inference(tree=None, seq_meta=None, field='country', confidence=True,
                        missing='?', sampling_bias_correction=None, weights=None):
    """
    Infer likely ancestral states of a discrete character assuming a time reversible model.

    Parameters
    ----------
    tree : str
        name of tree file
    seq_meta : dict
        meta data associated with sequences
    field : str, optional
        meta data field to use
    confidence : bool, optional
        calculate confidence values for inferences
    missing : str, optional
        character that is to be interpreted as missing data, default='?'
    sampling_bias_correction : None, optional
        factor by which the transition rate is scaled up to counter sampling bias
    weights : None, optional
        vector of equilibrium frequencies that one expects the far ancestor to be sampled from

    Returns
    -------
    T : Phylo.Tree
        Biophyton tree
    gtr : treetime.GTR
        GTR model
    alphabet : dict
        mapping of character states to
    """
    from treetime.wrappers import reconstruct_discrete_traits
    from Bio.Align import MultipleSeqAlignment
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio import Phylo

    T = Phylo.read(tree, 'newick')
    traits = {}
    nodes = {n.name:n for n in T.get_terminals()}
    for name, meta in seq_meta.items():
        if field in meta and name in nodes:
            traits[name] = meta[field]
    unique_states = list(set(traits.values()))

    if len(unique_states)==0:
        print("WARNING: no states found for discrete state reconstruction.")
        for node in T.find_clades():
            node.__setattr__(field, None)
        return T, None, {}
    elif len(unique_states)==1:
        print("WARNING: only one state found for discrete state reconstruction:", unique_states)
        for node in T.find_clades():
            node.__setattr__(field, unique_states[0])
        return T, None, {}
    elif len(unique_states)<300:
        tt, letter_to_state, reverse_alphabet = \
            reconstruct_discrete_traits(T, traits, missing_data=missing,
                 sampling_bias_correction=sampling_bias_correction, weights=weights)
    else:
        print("ERROR: 300 or more distinct discrete states found. TreeTime is currently not set up to handle that many states.")
        sys.exit(1)

    if tt is None:
        print("ERROR in discrete state reconstruction in TreeTime. Please look for errors above.")
        sys.exit(1)

    # attach inferred states as e.g. node.region = 'africa'
    for node in tt.tree.find_clades():
        node.__setattr__(field, letter_to_state[node.cseq[0]])

    # if desired, attach entropy and confidence as e.g. node.region_entropy = 0.03
    if confidence:
        for node in tt.tree.find_clades():
            pdis = node.marginal_profile[0]
            S = -np.sum(pdis*np.log(pdis+TINY))

            marginal = [(letter_to_state[tt.gtr.alphabet[i]], pdis[i]) for i in range(len(tt.gtr.alphabet))]
            marginal.sort(key=lambda x: x[1], reverse=True) # sort on likelihoods
            marginal = [(a, b) for a, b in marginal if b > 0.001][:4] #only take stuff over .1% and the top 4 elements
            conf = {a:b for a,b in marginal}
            node.__setattr__(field + "_entropy", S)
            node.__setattr__(field + "_confidence", conf)

    return tt.tree, tt.gtr, letter_to_state


def register_arguments(parser):
    """Add subcommand specific arguments

    Parameters
    ----------
    parser : argparse
        subcommand argument parser
    """
    parser.add_argument('--tree', '-t', required=True, help="tree to perform trait reconstruction on")
    parser.add_argument('--metadata', required=True, metavar="FILE", help="table with metadata, as CSV or TSV")
    parser.add_argument('--weights', required=False, help="tsv/csv table with equilibrium probabilities of discrete states")
    parser.add_argument('--columns', required=True, nargs='+',
                        help='metadata fields to perform discrete reconstruction on')
    parser.add_argument('--confidence',action="store_true",
                        help='record the distribution of subleading mugration states')
    parser.add_argument('--sampling-bias-correction', type=float,
                        help='a rough estimate of how many more events would have been observed'
                             ' if sequences represented an even sample. This should be'
                             ' roughly the (1-sum_i p_i^2)/(1-sum_i t_i^2), where p_i'
                             ' are the equilibrium frequencies and t_i are apparent ones.'
                             '(or rather the time spent in a particular state on the tree)')
    parser.add_argument('--output-node-data', type=str, help='name of JSON file to save trait inferences to')
    parser.epilog = "Note that missing data must be represented by a `?` character. Missing data will currently be inferred."

def run(args):
    """run mugration inference

    Parameters
    ----------
    args : namespace
        command line arguments are parsed by argparse
    """
    tree_fname = args.tree
    traits, columns = read_metadata(args.metadata)

    from Bio import Phylo
    T = Phylo.read(tree_fname, 'newick')
    missing_internal_node_names = [n.name is None for n in T.get_nonterminals()]
    if np.all(missing_internal_node_names):
        print("\n*** WARNING: Tree has no internal node names!")
        print("*** Without internal node names, ancestral traits can't be linked up to the correct node later.")
        print("*** If you want to use 'augur export' later, re-run this command with the output of 'augur refine'.")
        print("*** If you haven't run 'augur refine', you can add node names to your tree by running:")
        print("*** augur refine --tree %s --output-tree <filename>.nwk"%(tree_fname) )
        print("*** And use <filename>.nwk as the tree when running 'ancestral', 'translate', and 'traits'")

    if args.weights:
        weight_dict = {c:{} for c in args.columns}
        sep = ',' if args.weights.endswith('csv') else '\t'
        with open(args.weights, 'r', encoding='utf-8') as fh:
            for line in fh:
                if line[0]=='#':
                    continue
                name, trait, value = line.strip().split(sep)
                if name in weight_dict:
                    weight_dict[name][trait] = float(value)
        for c in weight_dict:
            if len(weight_dict[c])==0:
                weight_dict[c]=None
    else:
        weight_dict = {c:None for c in args.columns}

    mugration_states = defaultdict(dict)
    models = defaultdict(dict)
    out_prefix = '.'.join(args.output_node_data.split('.')[:-1])

    from treetime import version as treetime_version
    print(f"augur traits is using TreeTime version {treetime_version}")

    for column in args.columns:
        T, gtr, alphabet = mugration_inference(tree=tree_fname, seq_meta=traits,
                                               field=column, confidence=args.confidence,
                                               sampling_bias_correction=args.sampling_bias_correction,
                                               weights=weight_dict[column])
        if T is None: # something went wrong
            continue

        for node in T.find_clades():
            mugration_states[node.name][column] = node.__getattribute__(column)
            if args.confidence:
                mugration_states[node.name][column+'_confidence'] = node.__getattribute__(column+'_confidence')
                mugration_states[node.name][column+'_entropy'] = node.__getattribute__(column+'_entropy')

        if gtr:
            # add gtr models to json structure for export
            models[column]['rate'] = gtr.mu
            models[column]['alphabet'] = [alphabet[k] for k in sorted(alphabet.keys())]
            models[column]['equilibrium_probabilities'] = list(gtr.Pi)
            models[column]['transition_matrix'] = [list(x) for x in gtr.W]

        if gtr:
            with open(out_prefix+'%s.mugration_model.txt'%column, 'w', encoding='utf-8') as ofile:
                ofile.write('Map from character to field name\n')
                for k,v in alphabet.items():
                    ofile.write(k+':\t'+str(v)+'\n')
                ofile.write('\n\n')

                ofile.write(str(gtr))

    out_name = get_json_name(args, out_prefix+'_traits.json')
    write_json({"models":models, "nodes":mugration_states},out_name)

    print("\nInferred ancestral states of discrete character using TreeTime:"
          "\n\tSagulenko et al. TreeTime: Maximum-likelihood phylodynamic analysis"
          "\n\tVirus Evolution, vol 4, https://academic.oup.com/ve/article/4/1/vex042/4794731\n", file=sys.stdout)

    print("results written to", out_name, file=sys.stdout)
