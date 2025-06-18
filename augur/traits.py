"""
Infer ancestral traits based on a tree.
"""

import numpy as np
from collections import defaultdict, OrderedDict, Counter
import sys
from .argparse_ import ExtendOverwriteDefault
from .errors import AugurError
from .io.file import open_file
from .io.metadata import DEFAULT_DELIMITERS, DEFAULT_ID_COLUMNS, InvalidDelimiter, read_metadata
from .utils import write_json, get_json_name
TINY = 1e-12

def mugration_inference(tree=None, seq_meta=None, field='country', confidence=True,
                        missing='?', sampling_bias_correction=None, weights=None):
    """
    Infer likely ancestral states of a discrete character assuming a time reversible model.

    Parameters
    ----------
    tree : str
        name of tree file
    seq_meta : pandas.DataFrame
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
    T : Bio.Phylo.BaseTree.Tree
        Biophyton tree
    gtr : treetime.GTR
        GTR model
    alphabet : dict
        mapping of character states to
    """
    from treetime.wrappers import reconstruct_discrete_traits
    from Bio import Phylo

    T = Phylo.read(tree, 'newick')
    traits = {}
    nodes = {n.name:n for n in T.get_terminals()}
    for name, meta in seq_meta.iterrows():
        if field in meta and name in nodes and meta[field] != missing:
            traits[name] = meta[field]
    unique_states = list(set(traits.values()))

    if len(unique_states)==0:
        print("WARNING: no states found for discrete state reconstruction.", file=sys.stderr)
        for node in T.find_clades():
            node.__setattr__(field, None)
        return T, None, {}
    elif len(unique_states)==1:
        print("WARNING: only one state found for discrete state reconstruction:", unique_states, file=sys.stderr)
        for node in T.find_clades():
            node.__setattr__(field, unique_states[0])
        return T, None, {}
    elif len(unique_states)<300:
        tt, letter_to_state, reverse_alphabet = \
            reconstruct_discrete_traits(T, traits, missing_data=missing,
                 sampling_bias_correction=sampling_bias_correction, weights=weights)
    else:
        print("ERROR: 300 or more distinct discrete states found. TreeTime is currently not set up to handle that many states.", file=sys.stderr)
        sys.exit(1)

    if tt is None:
        print("ERROR in discrete state reconstruction in TreeTime. Please look for errors above.", file=sys.stderr)
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
            # Values are defined for all demes, although many/most will be 0, so only take those over 0.1%
            marginal = [(a, b) for a, b in marginal if b > 0.001]
            conf = {a:b for a,b in marginal}
            node.__setattr__(field + "_entropy", S)
            node.__setattr__(field + "_confidence", conf)

    return tt.tree, tt.gtr, letter_to_state


class BranchLabeller():
    """
    A class to create branch labels based on changes in *column* state.
    If the *enabled* arg is false then the user-facing methods are no-ops.
    """
    def __init__(self, labels_arg, confidence_arg, column):
        self.labels = {}

        # parse the label argument
        self.column = None
        for word in labels_arg:
            parts = word.split("=")
            if len(parts)>2:
                raise AugurError("The --branch-labels argument {labels_arg:r} is not valid (multiple equals signs in a value)") 
            if parts[0]==column:
                self.column = column
                self.column_label = parts[1] if len(parts)==2 else column

        if self.column is None:
            return

        # parse the confidence argument
        self._confidence_threshold = None
        for word in confidence_arg:
            parts = word.split("=")
            if len(parts)!=2:
                raise AugurError("The --branch-confidence argument {confidence_arg:r} is not valid (each element must be COLUMN=CONFIDENCE)") 
            if parts[0]==column:
                self._parse_confidence(parts[1])

    def _parse_confidence(self, value): # simpler than using @property
        if len(value): # an empty value indicates no confidence threshold
            try:
                self._confidence_threshold = int(value)
                assert self._confidence_threshold >=0
                assert self._confidence_threshold <=100
            except:
                raise AugurError(f"--branch-labels' confidence threshold must be an integer between 0 and 100  (column {self.column})")

    def _state(self, node):
        attr = getattr(node, self.column)
        if self._confidence_threshold is None:
            return attr
        conf = getattr(node, f"{self.column}_confidence", None)[attr]*100
        if conf > self._confidence_threshold:
            return attr
        return "uncertain"

    def process(self, node):
        if self.column is None:
            return
        node_state = self._state(node)
        if not node.up:
            self.labels[node.name] = node_state
            return
        if (parent_state:=self._state(node.up)) != node_state:
            self.labels[node.name] = f"{parent_state} → {node_state}"
        return

    def changes(self):
        if self.column is None:
            return
        counts = Counter(self.labels.values())
        observed = defaultdict(int)
        for node_name, base_label in self.labels.items():
            if counts[base_label]==1:
                label = base_label
            else:
                observed[base_label]+=1
                label = f"{base_label} {observed[base_label]}"
            yield (self.column_label, node_name, label)

def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("traits", help=__doc__)
    parser.add_argument('--tree', '-t', required=True, help="tree to perform trait reconstruction on")
    parser.add_argument('--metadata', required=True, metavar="FILE", help="table with metadata")
    parser.add_argument('--metadata-delimiters', default=DEFAULT_DELIMITERS, nargs="+", action=ExtendOverwriteDefault,
                        help="delimiters to accept when reading a metadata file. Only one delimiter will be inferred.")
    parser.add_argument('--metadata-id-columns', default=DEFAULT_ID_COLUMNS, nargs="+", action=ExtendOverwriteDefault,
                        help="names of possible metadata columns containing identifier information, ordered by priority. Only one ID column will be inferred.")
    parser.add_argument('--weights', required=False, help="tsv/csv table with equilibrium probabilities of discrete states")
    parser.add_argument('--columns', required=True, nargs='+', action=ExtendOverwriteDefault,
                        help='metadata fields to perform discrete reconstruction on')
    parser.add_argument('--confidence',action="store_true",
                        help='record the distribution of subleading mugration states')
    parser.add_argument('--branch-labels', nargs='+', default=[], metavar="COLUMN[=NAME]", action=ExtendOverwriteDefault,
                        help='Add branch labels where there is a change in trait inferred for that column.' \
                        ' You must supply this for each column you would like to label.' \
                        ' By default the branch label key the same as the column name, but you may customise this' \
                        ' via the COLUMN=NAME syntax.')
    parser.add_argument('--branch-confidence', nargs='+', default=[], metavar="COLUMN=CONFIDENCE", action=ExtendOverwriteDefault,
                        help='Only label state changes where the confidence percentage is above the specified value.' \
                        'Transitions to lower confidence states will be represented by a "uncertain" label.')
    parser.add_argument('--sampling-bias-correction', type=float,
                        help='a rough estimate of how many more events would have been observed'
                             ' if sequences represented an even sample. This should be'
                             ' roughly the (1-sum_i p_i^2)/(1-sum_i t_i^2), where p_i'
                             ' are the equilibrium frequencies and t_i are apparent ones.'
                             '(or rather the time spent in a particular state on the tree)')
    parser.add_argument('--output-node-data', type=str, help='name of JSON file to save trait inferences to')
    parser.epilog = "Note that missing data must be represented by a `?` character. Missing data will currently be inferred."
    return parser


def run(args):
    """run mugration inference

    Parameters
    ----------
    args : argparse.Namespace
        command line arguments are parsed by argparse
    """
    tree_fname = args.tree
    try:
        traits = read_metadata(
            args.metadata,
            delimiters=args.metadata_delimiters,
            id_columns=args.metadata_id_columns,

            # Read all columns as string for discrete trait analysis
            dtype="string",
        )
    except InvalidDelimiter:
        raise AugurError(
                f"Could not determine the delimiter of {args.metadata!r}. "
                f"Valid delimiters are: {args.metadata_delimiters!r}. "
                "This can be changed with --metadata-delimiters."
            )

    from Bio import Phylo
    T = Phylo.read(tree_fname, 'newick')
    missing_internal_node_names = [n.name is None for n in T.get_nonterminals()]
    if np.all(missing_internal_node_names):
        print("\n*** WARNING: Tree has no internal node names!", file=sys.stderr)
        print("*** Without internal node names, ancestral traits can't be linked up to the correct node later.", file=sys.stderr)
        print("*** If you want to use 'augur export' later, re-run this command with the output of 'augur refine'.", file=sys.stderr)
        print("*** If you haven't run 'augur refine', you can add node names to your tree by running:", file=sys.stderr)
        print("*** augur refine --tree %s --output-tree <filename>.nwk" % (tree_fname), file=sys.stderr)
        print("*** And use <filename>.nwk as the tree when running 'ancestral', 'translate', and 'traits'", file=sys.stderr)

    if args.weights:
        weight_dict = {c:{} for c in args.columns}
        sep = ',' if args.weights.endswith('csv') else '\t'
        with open_file(args.weights, 'r') as fh:
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
    branch_states = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
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

        branch_labeller = BranchLabeller(args.branch_labels, args.branch_confidence, column)

        for node in T.find_clades():
            mugration_states[node.name][column] = getattr(node, column)

            if args.confidence:
                confidence = getattr(node, f"{column}_confidence", None)
                if confidence is not None:
                    mugration_states[node.name][f"{column}_confidence"] = confidence

                entropy = getattr(node, f"{column}_entropy", None)
                if entropy is not None:
                    mugration_states[node.name][f"{column}_entropy"] = entropy

            branch_labeller.process(node)

        for column_name, node_name, label in branch_labeller.changes():
            branch_states[node_name]['labels'][column_name] = label

        if gtr:
            # add gtr models to json structure for export
            models[column]['rate'] = gtr.mu
            models[column]['alphabet'] = [alphabet[k] for k in sorted(alphabet.keys())]
            models[column]['equilibrium_probabilities'] = list(gtr.Pi)
            models[column]['transition_matrix'] = [list(x) for x in gtr.W]

        if gtr:
            with open_file(out_prefix+'%s.mugration_model.txt'%column, 'w') as ofile:
                ofile.write('Map from character to field name\n')
                for k,v in alphabet.items():
                    ofile.write(k+':\t'+str(v)+'\n')
                ofile.write('\n\n')

                ofile.write(str(gtr))

    out_name = get_json_name(args, out_prefix+'_traits.json')
    json_data = OrderedDict([["models", models], ["nodes", mugration_states]])
    if branch_states:
        json_data['branches'] = branch_states
    write_json(json_data, out_name)

    print("\nInferred ancestral states of discrete character using TreeTime:"
          "\n\tSagulenko et al. TreeTime: Maximum-likelihood phylodynamic analysis"
          "\n\tVirus Evolution, vol 4, https://academic.oup.com/ve/article/4/1/vex042/4794731\n", file=sys.stdout)

    print("results written to", out_name, file=sys.stdout)
