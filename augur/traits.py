import numpy as np
from collections import defaultdict
import json, os
import pandas as pd
from .utils import read_metadata
TINY = 1e-12

def mugration_inference(tree=None, seq_meta=None, field='country', confidence=True,
                        infer_gtr=True, root_state=None, missing='?'):
    from treetime import GTR
    from Bio.Align import MultipleSeqAlignment
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio import Phylo

    T = Phylo.read(tree, 'newick')
    nodes = {n.name:n for n in T.get_terminals()}

    # Determine alphabet only counting tips in the tree
    places = set()
    for name, meta in seq_meta.items():
        if field in meta and name in nodes:
            places.add(meta[field])
    if root_state is not None:
        places.add(root_state)

    # construct GTR (flat for now). The missing DATA symbol is a '-' (ord('-')==45)
    places = sorted(places)
    nc = len(places)
    if nc>180:
        print("ERROR: geo_inference: can't have more than 180 places!")
        return None,None
    elif nc==1:
        print("WARNING: geo_inference: only one place found -- set every internal node to %s!"%places[0])
        return None,None
    elif nc==0:
        print("ERROR: geo_inference: list of places is empty!")
        return None,None
    else:
        # set up model
        alphabet = {chr(65+i):place for i,place in enumerate(places)}
        model = GTR.custom(pi = np.ones(nc, dtype=float)/nc, W=np.ones((nc,nc)),
                          alphabet = np.array(sorted(alphabet.keys())))

        missing_char = chr(65+nc)
        alphabet[missing_char]=missing
        model.profile_map[missing_char] = np.ones(nc)
        model.ambiguous = missing_char
        alphabet_rev = {v:k for k,v in alphabet.items()}

        # construct pseudo alignment
        pseudo_seqs = []
        for name, meta in seq_meta.items():
            if name in nodes:
                s=alphabet_rev[meta[field]] if field in meta else missing_char
                pseudo_seqs.append(SeqRecord(Seq(s), name=name, id=name))
        aln = MultipleSeqAlignment(pseudo_seqs)

        # set up treetime and infer
        from treetime import TreeAnc
        tt = TreeAnc(tree=tree, aln=aln, gtr=model, convert_upper=False, verbose=0)
        tt.use_mutation_length=False
        tt.infer_ancestral_sequences(infer_gtr=infer_gtr, store_compressed=False, pc=5.0,
                                     marginal=True, normalized_rate=False)

        # attach inferred states as e.g. node.region = 'africa'
        for node in tt.tree.find_clades():
            node.__setattr__(field, alphabet[node.sequence[0]])

        # if desired, attach entropy and confidence as e.g. node.region_entropy = 0.03
        if confidence:
            for node in tt.tree.find_clades():
                pdis = node.marginal_profile[0]
                S = -np.sum(pdis*np.log(pdis+TINY))

                marginal = [(alphabet[tt.gtr.alphabet[i]], pdis[i]) for i in range(len(tt.gtr.alphabet))]
                marginal.sort(key=lambda x: x[1], reverse=True) # sort on likelihoods
                marginal = [(a, b) for a, b in marginal if b > 0.001][:4] #only take stuff over .1% and the top 4 elements
                conf = {a:b for a,b in marginal}
                node.__setattr__(field + "_entropy", S)
                node.__setattr__(field + "_confidence", conf)

        return tt, alphabet



def run(args):
    tree_fname = args.tree
    traits, columns = read_metadata(args.metadata)

    mugration_states = defaultdict(dict)
    for column in args.columns:
        tt, alphabet = mugration_inference(tree=tree_fname, seq_meta=traits,
                            field=column, confidence=args.confidence)
        if tt is None: # something went wrong
            continue


        for node in tt.tree.find_clades():
            mugration_states[node.name][column] = node.__getattribute__(column)
            if args.confidence:
                mugration_states[node.name][column+'_confidence'] = node.__getattribute__(column+'_confidence')
                mugration_states[node.name][column+'_entropy'] = node.__getattribute__(column+'_entropy')

        #if args.output is default (no dir), including '/' messes up writing
        prefix = os.path.dirname(args.output)+'/' if len(os.path.dirname(args.output)) != 0 else ''
        with open(prefix+'%s.mugration_model.txt'%column, 'w') as ofile:
            ofile.write('Map from character to field name\n')
            for k,v in alphabet.items():
                ofile.write(k+':\t'+str(v)+'\n')
            ofile.write('\n\n')

            ofile.write(str(tt.gtr))

    with open(args.output, 'w') as results:
        json.dump({"nodes":mugration_states}, results, indent=1)

    print("\nInferred ancestral states of discrete character using TreeTime:"
          "\n\tSagulenko et al. TreeTime: Maximum-likelihood phylodynamic analysis"
          "\n\tVirus Evolution, vol 4, https://academic.oup.com/ve/article/4/1/vex042/4794731\n")
    print("results written to",args.output)
