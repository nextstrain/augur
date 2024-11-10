"""
Refine an initial tree using sequence metadata.
"""
import numpy as np
import sys
from Bio import Phylo
from .argparse_ import ExtendOverwriteDefault
from .dates import get_numerical_dates
from .dates.errors import InvalidYearBounds
from .io.metadata import DEFAULT_DELIMITERS, DEFAULT_ID_COLUMNS, METADATA_DATE_COLUMN, InvalidDelimiter, Metadata, read_metadata
from .utils import read_tree, write_json, InvalidTreeError
from .errors import AugurError
from treetime.vcf_utils import read_vcf
from treetime.seq_utils import profile_maps

def refine(tree=None, aln=None, ref=None, dates=None, branch_length_inference='auto',
             confidence=False, resolve_polytomies=True, stochastic_resolve=False, max_iter=2, precision='auto',
             infer_gtr=True, Tc=0.01, reroot=None, use_marginal='always', fixed_pi=None, use_fft=True,
             clock_rate=None, clock_std=None, clock_filter_iqd=None, verbosity=1, covariance=True, rng_seed=None, **kwarks):
    from treetime import TreeTime

    try: #Tc could be a number or  'opt' or 'skyline'. TreeTime expects a float or int if a number.
        Tc = float(Tc)
    except ValueError:
        True #let it remain a string

    if (ref is not None) and (fixed_pi is None): #if VCF, fix pi
        #Otherwise mutation TO gaps is overestimated b/c of seq length
        fixed_pi = [ref.count(base)/len(ref) for base in ['A','C','G','T','-']]
        if fixed_pi[-1] == 0:
            fixed_pi[-1] = 0.05
            fixed_pi = [v-0.01 for v in fixed_pi]

    if ref is not None: # VCF -> adjust branch length
        #set branch length mode explicitly if auto, as informative-site only
        #trees can have big branch lengths, making this set incorrectly in TreeTime
        if branch_length_inference == 'auto':
            branch_length_inference = 'joint'

    #send ref, if is None, does no harm
    tt = TreeTime(tree=tree, aln=aln, ref=ref, dates=dates, use_fft=use_fft,
                  verbose=verbosity, gtr='JC69', precision=precision, rng_seed=rng_seed)

    # conditionally run clock-filter and remove bad tips
    if clock_filter_iqd:
        # treetime clock filter will mark, but not remove bad tips
        tt.clock_filter(reroot=reroot, n_iqd=clock_filter_iqd, plot=False) #use whatever was specified
        # remove them explicitly
        leaves = [x for x in tt.tree.get_terminals()]
        for n in leaves:
            if n.bad_branch:
                tt.tree.prune(n)
                print('pruning leaf ', n.name)
        # fix treetime set-up for new tree topology
        tt.prepare_tree()

    if confidence and use_marginal:
        # estimate confidence intervals via marginal ML and assign
        # marginal ML times to nodes
        marginal = 'always' if use_fft else 'assign'
    else:
        marginal = 'only-final' if confidence else False

    # uncertainty of the the clock rate is relevant if confidence intervals are estimated
    if confidence and clock_std:
        vary_rate = clock_std # if standard devivation of clock is specified, use that
    elif (clock_rate is None) and confidence and covariance:
        vary_rate = True      # if run in covariance mode, standard deviation can be estimated
    else:
        vary_rate = False     # otherwise, rate uncertainty will be ignored

    tt.run(infer_gtr=infer_gtr, root=reroot, Tc=Tc, time_marginal=marginal,
           branch_length_mode=branch_length_inference, resolve_polytomies=resolve_polytomies,
           stochastic_resolve=stochastic_resolve, max_iter=max_iter, fixed_pi=fixed_pi,
           fixed_clock_rate=clock_rate, vary_rate=vary_rate, use_covariation=covariance,
           raise_uncaught_exceptions=True, **kwarks)

    if confidence:
        for n in tt.tree.find_clades():
            n.num_date_confidence = list(tt.get_max_posterior_region(n, 0.9))

    print("\nInferred a time resolved phylogeny using TreeTime:"
          "\n\tSagulenko et al. TreeTime: Maximum-likelihood phylodynamic analysis"
          "\n\tVirus Evolution, vol 4, https://academic.oup.com/ve/article/4/1/vex042/4794731\n")
    return tt


def collect_node_data(T, attributes):
    data = {}
    for n in T.find_clades():
        data[n.name] = {
            attr: getattr(n, attr)
            for attr in attributes
            if getattr(n, attr, None) is not None
        }
    return data


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("refine", help=__doc__)
    parser.add_argument('--alignment', '-a', help="alignment in fasta or VCF format")
    parser.add_argument('--tree', '-t', required=True, help="prebuilt Newick")
    parser.add_argument('--metadata', type=str, metavar="FILE", help="sequence metadata")
    parser.add_argument('--metadata-delimiters', default=DEFAULT_DELIMITERS, nargs="+", action=ExtendOverwriteDefault,
                        help="delimiters to accept when reading a metadata file. Only one delimiter will be inferred.")
    parser.add_argument('--metadata-id-columns', default=DEFAULT_ID_COLUMNS, nargs="+", action=ExtendOverwriteDefault,
                        help="names of possible metadata columns containing identifier information, ordered by priority. Only one ID column will be inferred.")
    parser.add_argument('--output-tree', type=str, help='file name to write tree to')
    parser.add_argument('--output-node-data', type=str, help='file name to write branch lengths as node data')
    parser.add_argument('--use-fft', action="store_true", help="produce timetree using FFT for convolutions")
    parser.add_argument('--max-iter', default=2, type=int, help="maximal number of iterations TreeTime uses for timetree inference")
    parser.add_argument('--timetree', action="store_true", help="produce timetree using treetime, requires tree where branch length is in units of average number of nucleotide or protein substitutions per site (and branch lengths do not exceed 4)")
    parser.add_argument('--coalescent', help="coalescent time scale in units of inverse clock rate (float), optimize as scalar ('opt'), or skyline ('skyline')")
    parser.add_argument('--gen-per-year', default=50, type=float, help="number of generations per year, relevant for skyline output('skyline')")
    parser.add_argument('--clock-rate', type=float, help="fixed clock rate")
    parser.add_argument('--clock-std-dev', type=float, help="standard deviation of the fixed clock_rate estimate")
    parser.add_argument('--root', nargs="+", action=ExtendOverwriteDefault, default='best', help="rooting mechanism ('best', least-squares', 'min_dev', 'oldest', 'mid_point') "
                                "OR node to root by OR two nodes indicating a monophyletic group to root by. "
                                "Run treetime -h for definitions of rooting methods.")
    parser.add_argument('--keep-root', action="store_true", help="do not reroot the tree; use it as-is. "
                                "Overrides anything specified by --root.")
    parser.add_argument('--covariance', dest='covariance', action='store_true', help="Account for covariation when estimating "
                                "rates and/or rerooting. "
                                "Use --no-covariance to turn off.")
    parser.add_argument('--no-covariance', dest='covariance', action='store_false')  #If you set help here, it displays 'default: True' - which is confusing!

    resolve_group = parser.add_mutually_exclusive_group()
    resolve_group.add_argument('--keep-polytomies', action='store_true', help='Do not attempt to resolve polytomies')
    resolve_group.add_argument('--stochastic-resolve', action='store_true', help='Resolve polytomies via stochastic subtree building rather than greedy optimization')
    resolve_group.add_argument('--greedy-resolve', action='store_false', dest='stochastic_resolve') # inverse of `--stochastic-resolve` to facilitate changing defaults in the future

    parser.add_argument('--precision', type=int, choices=[0,1,2,3], help="precision used by TreeTime to determine the number of grid points that are used for the evaluation of the branch length interpolation objects. Values range from 0 (rough) to 3 (ultra fine) and default to 'auto'.")
    parser.add_argument('--date-format', default="%Y-%m-%d", help="date format")
    parser.add_argument('--date-confidence', action="store_true", help="calculate confidence intervals for node dates")
    parser.add_argument('--date-inference', default='joint', choices=["joint", "marginal"],
                                help="assign internal nodes to their marginally most likely dates, not jointly most likely")
    parser.add_argument('--branch-length-inference', default='auto', choices = ['auto', 'joint', 'marginal', 'input'],
                                help='branch length mode of treetime to use')
    parser.add_argument('--clock-filter-iqd', type=float, help='clock-filter: remove tips that deviate more than n_iqd '
                                'interquartile ranges from the root-to-tip vs time regression')
    parser.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the VCF was mapped to')
    parser.add_argument('--year-bounds', type=int, nargs='+', action='extend', help='specify min or max & min prediction bounds for samples with XX in year')
    parser.add_argument('--divergence-units', type=str, choices=['mutations', 'mutations-per-site'],
                        default='mutations-per-site', help='Units in which sequence divergences is exported.')
    parser.add_argument('--seed', type=int, help='seed for random number generation')
    parser.add_argument('--verbosity', type=int, default=1, help='treetime verbosity, between 0 and 6 (higher values more output)')
    parser.set_defaults(covariance=True)
    return parser


def run(args):

    # check alignment type, set flags, read in if VCF
    is_vcf = False
    ref = None

    # node data is the dict that will be exported as json
    node_data = {'alignment': args.alignment}
    # list of node attributes that are to be exported, will grow
    attributes = ['branch_length', 'confidence']

    try:
        T = read_tree(args.tree)
        node_data['input_tree'] = args.tree
    except (FileNotFoundError, InvalidTreeError) as error:
        print("ERROR: %s" % error, file=sys.stderr)
        return 1

    if not args.alignment:
        if args.timetree:
            print("ERROR: alignment is required for ancestral reconstruction or timetree inference", file=sys.stderr)
            return 1

        if args.divergence_units=='mutations':
            print("ERROR: alignment is required for divergence in units of mutations", file=sys.stderr)
            return 1

        # fake alignment to appease treetime when only using it for naming nodes...
        from Bio import SeqRecord, Seq, Align
        seqs = []
        for n in T.get_terminals():
            seqs.append(SeqRecord.SeqRecord(seq=Seq.Seq('ACGT'), id=n.name, name=n.name, description=''))
        aln = Align.MultipleSeqAlignment(seqs)
    elif any([args.alignment.lower().endswith(x) for x in ['.vcf', '.vcf.gz']]):
        if not args.vcf_reference:
            print("ERROR: a reference Fasta is required with VCF-format alignments", file=sys.stderr)
            return 1

        compress_seq = read_vcf(args.alignment, args.vcf_reference)
        aln = compress_seq['sequences']
        ref = compress_seq['reference']
        is_vcf = True
    else:
        aln = args.alignment

    from treetime import version as treetime_version
    print(f"augur refine is using TreeTime version {treetime_version}")

    # if not specified, construct default output file name with suffix _tt.nwk
    if args.output_tree:
        tree_fname = args.output_tree
    elif args.alignment:
        tree_fname = '.'.join(args.alignment.split('.')[:-1]) + '_tt.nwk'
    else:
        tree_fname = '.'.join(args.tree.split('.')[:-1]) + '_tt.nwk'

    if args.root and len(args.root) == 1: #if anything but a list of seqs, don't send as a list
        args.root = args.root[0]
    if args.keep_root:  # This flag overrides anything specified by 'root'
        args.root = None

    if args.timetree:
        # load meta data and covert dates to numeric
        if args.metadata is None:
            print("ERROR: meta data with dates is required for time tree reconstruction", file=sys.stderr)
            return 1

        try:
            metadata_object = Metadata(args.metadata, args.metadata_delimiters, args.metadata_id_columns)
        except InvalidDelimiter:
            raise AugurError(
                f"Could not determine the delimiter of {args.metadata!r}. "
                f"Valid delimiters are: {args.metadata_delimiters!r}. "
                "This can be changed with --metadata-delimiters."
            )

        metadata = read_metadata(
            args.metadata,
            delimiters=[metadata_object.delimiter],
            columns=[metadata_object.id_column, METADATA_DATE_COLUMN],
            id_columns=[metadata_object.id_column],
            dtype="string",
        )

        try:
            dates = get_numerical_dates(metadata, fmt=args.date_format,
                                        min_max_year=args.year_bounds)
        except InvalidYearBounds as error:
            raise AugurError(f"Invalid value for --year-bounds: {error}")

        # save input state string for later export
        for n in T.get_terminals():
            if n.name in metadata.index and METADATA_DATE_COLUMN in metadata.columns:
                n.raw_date = metadata.at[n.name, METADATA_DATE_COLUMN]

        if args.date_confidence:
            time_inference_mode = 'always' if args.date_inference=='marginal' else 'only-final'
        else:
            time_inference_mode = 'always' if args.date_inference=='marginal' else 'never'

        if args.root == 'mid_point':
            # root at midpoint and disable downstream rerooting in TreeTime
            T.root_at_midpoint()
            args.root = None

        tt = refine(tree=T, aln=aln, ref=ref, dates=dates, confidence=args.date_confidence,
                    reroot=args.root, # or 'best', # We now have a default in param spec - this just adds confusion.
                    Tc=0.01 if args.coalescent is None else args.coalescent, #use 0.01 as default coalescent time scale
                    use_marginal = time_inference_mode, use_fft=args.use_fft,
                    branch_length_inference = args.branch_length_inference or 'auto',
                    precision = 'auto' if args.precision is None else args.precision,
                    clock_rate=args.clock_rate, clock_std=args.clock_std_dev,
                    clock_filter_iqd=args.clock_filter_iqd, max_iter=args.max_iter,
                    covariance=args.covariance, resolve_polytomies=(not args.keep_polytomies),
                    stochastic_resolve=args.stochastic_resolve, verbosity=args.verbosity, rng_seed=args.seed)

        node_data['clock'] = {'rate': tt.date2dist.clock_rate,
                              'intercept': tt.date2dist.intercept,
                              'rtt_Tmrca': -tt.date2dist.intercept/tt.date2dist.clock_rate}
        # Include the standard deviation of the clock rate, if the covariance
        # matrix is available.
        if hasattr(tt.date2dist, "cov") and tt.date2dist.cov is not None:
            node_data["clock"]["cov"] = tt.date2dist.cov
            node_data["clock"]["rate_std"] = np.sqrt(tt.date2dist.cov[0, 0])

        if args.coalescent=='skyline':
            try:
                skyline, conf = tt.merger_model.skyline_inferred(gen=args.gen_per_year, confidence=2)
                node_data['skyline'] = [[float(x) for x in skyline.x], [float(y) for y in conf[0]],
                                        [float(y) for y in skyline.y], [float(y) for y in conf[1]]]
            except:
                print(
                    "ERROR: skyline optimization by TreeTime has failed.",
                    "To avoid this error, try running without coalescent optimization or with `--coalescent opt` instead of `--coalescent skyline`.",
                    file=sys.stderr
                )
                return 1

        attributes.extend(['numdate', 'clock_length', 'mutation_length', 'raw_date', 'date'])
        if args.date_confidence:
            attributes.append('num_date_confidence')
    else:
        from treetime import TreeAnc
        # instantiate treetime for the sole reason to name internal nodes
        if args.root:
            if args.root == 'best':
                print("Warning: To root without inferring a timetree, you must specify an explicit outgroup.")
                print("\tProceeding without re-rooting. To suppress this message, use '--keep-root'.\n")
            elif args.root in ['least-squares', 'oldest', 'min_dev']:
                raise TypeError("The rooting option '%s' is only available when inferring a timetree. Please specify an explicit outgroup."%args.root)
            elif args.root=="mid_point":
                T.root_at_midpoint()
            else:
                try:
                    T.root_with_outgroup(args.root)
                except ValueError as err:
                    raise ValueError(f"HINT: This error may be because your specified root with name '{args.root}' was not found in your alignment file") from err

        tt = TreeAnc(tree=T, aln=aln, ref=ref, gtr='JC69', verbose=args.verbosity, rng_seed=args.seed)

    node_data['nodes'] = collect_node_data(T, attributes)
    if args.divergence_units=='mutations-per-site': #default
        pass
    elif args.divergence_units=='mutations':
        if not args.timetree:
            tt.infer_ancestral_sequences()
        nuc_map = profile_maps['nuc']

        def are_sequence_states_different(nuc1, nuc2):
            '''
            determine whether two ancestral states should count as mutation for divergence estimates
            while correctly accounting for ambiguous nucleotides
            '''
            if nuc1 in ['-', 'N'] or nuc2 in ['-', 'N']:
                return False
            elif nuc1 in nuc_map and nuc2 in nuc_map:
                return np.sum(nuc_map[nuc1]*nuc_map[nuc2])==0
            else:
                return False

        for node in T.find_clades():
            n_muts = len([
                position
                for ancestral, position, derived in node.mutations
                if are_sequence_states_different(ancestral, derived)
            ])

            if args.timetree:
                node_data['nodes'][node.name]['mutation_length'] = n_muts

            node_data['nodes'][node.name]['branch_length'] = n_muts
    else:
        print("ERROR: divergence unit",args.divergence_units,"not supported!", file=sys.stderr)
        return 1

    # Export refined tree and node data
    tree_success = Phylo.write(T, tree_fname, 'newick', format_branch_length='%1.8f', branch_length_only=True)
    print("updated tree written to",tree_fname, file=sys.stdout)

    if args.output_node_data:
        node_data_fname = args.output_node_data
    elif args.alignment:
        node_data_fname = '.'.join(args.alignment.split('.')[:-1]) + '.node_data.json'
    else:
        node_data_fname = '.'.join(args.tree.split('.')[:-1]) + '.node_data.json'

    write_json(node_data, node_data_fname)
    print("node attributes written to",node_data_fname, file=sys.stdout)

    return 0 if tree_success else 1
