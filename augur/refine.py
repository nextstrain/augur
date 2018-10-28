"""
Refine an initial tree using sequence metadata.
"""

import os, shutil, time, sys
from Bio import Phylo
from .utils import read_metadata, get_numerical_dates, write_json
from treetime.vcf_utils import read_vcf, write_vcf

def refine(tree=None, aln=None, ref=None, dates=None, branch_length_inference='auto',
             confidence=False, resolve_polytomies=True, max_iter=2,
             infer_gtr=True, Tc=0.01, reroot=None, use_marginal=False, fixed_pi=None,
             clock_rate=None, clock_std=None, clock_filter_iqd=None, verbosity=1, **kwarks):
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
    tt = TreeTime(tree=tree, aln=aln, ref=ref, dates=dates,
                  verbose=verbosity, gtr='JC69')

    # conditionally run clock-filter and remove bad tips
    if clock_filter_iqd:
        # treetime clock filter will mark, but not remove bad tips
        tt.clock_filter(reroot='best', n_iqd=clock_filter_iqd, plot=False)
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
        marginal = 'assign'
    else:
        marginal = confidence

    vary_rate = False
    if clock_rate and clock_std:
        vary_rate = clock_std
    else:
        vary_rate = True

    tt.run(infer_gtr=infer_gtr, root=reroot, Tc=Tc, time_marginal=marginal,
           branch_length_mode=branch_length_inference, resolve_polytomies=resolve_polytomies,
           max_iter=max_iter, fixed_pi=fixed_pi, fixed_clock_rate=clock_rate,
           vary_rate=vary_rate, **kwarks)

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
        data[n.name] = {attr:n.__getattribute__(attr)
                        for attr in attributes if hasattr(n,attr)}
    return data


def register_arguments(parser):
    parser.add_argument('--alignment', '-a', help="alignment in fasta or VCF format")
    parser.add_argument('--tree', '-t', required=True, help="prebuilt Newick")
    parser.add_argument('--metadata', type=str, help="tsv/csv table with meta data for sequences")
    parser.add_argument('--output-tree', type=str, help='file name to write tree to')
    parser.add_argument('--output-node-data', type=str, help='file name to write branch lengths as node data')
    parser.add_argument('--timetree', action="store_true", help="produce timetree using treetime")
    parser.add_argument('--coalescent', help="coalescent time scale in units of inverse clock rate (float), optimize as scalar ('opt'), or skyline ('skyline')")
    parser.add_argument('--clock-rate', type=float, help="fixed clock rate")
    parser.add_argument('--clock-std-dev', type=float, help="standard deviation of the fixed clock_rate estimate")
    parser.add_argument('--root', nargs="+", help="rooting mechanism ('best', 'residual', 'rsq', 'min_dev') "
                                "OR node to root by OR two nodes indicating a monophyletic group to root by")
    parser.add_argument('--date-format', default="%Y-%m-%d", help="date format")
    parser.add_argument('--date-confidence', action="store_true", help="calculate confidence intervals for node dates")
    parser.add_argument('--date-inference', default='joint', choices=["joint", "marginal"],
                                help="assign internal nodes to their marginally most likely dates, not jointly most likely")
    parser.add_argument('--branch-length-inference', default='auto', choices = ['auto', 'joint', 'marginal', 'input'],
                                help='branch length mode of treetime to use')
    parser.add_argument('--clock-filter-iqd', type=float, help='clock-filter: remove tips that deviate more than n_iqd '
                                'interquartile ranges from the root-to-tip vs time regression')
    parser.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the VCF was mapped to')
    parser.add_argument('--year-bounds', type=int, nargs='+', help='specify min or max & min prediction bounds for samples with XX in year')


def run(args):
    # check alignment type, set flags, read in if VCF
    is_vcf = False
    ref = None

    # node data is the dict that will be exported as json
    node_data = {'alignment': args.alignment}
    # list of node attributes that are to be exported, will grow
    attributes = ['branch_length']

    # check if tree is provided an can be read
    for fmt in ["newick", "nexus"]:
        try:
            T = Phylo.read(args.tree, fmt)
            node_data['input_tree'] = args.tree
            break
        except:
            pass
    if T is None:
        print("ERROR: reading tree from %s failed."%args.tree)
        return 1

    if not args.alignment:
        # fake alignment to appease treetime when only using it for naming nodes...
        if args.timetree:
            print("ERROR: alignment is required for ancestral reconstruction or timetree inference")
            return 1
        from Bio import SeqRecord, Seq, Align
        seqs = []
        for n in T.get_terminals():
            seqs.append(SeqRecord.SeqRecord(seq=Seq.Seq('ACGT'), id=n.name, name=n.name, description=''))
        aln = Align.MultipleSeqAlignment(seqs)
    elif any([args.alignment.lower().endswith(x) for x in ['.vcf', '.vcf.gz']]):
        if not args.vcf_reference:
            print("ERROR: a reference Fasta is required with VCF-format alignments")
            return 1

        compress_seq = read_vcf(args.alignment, args.vcf_reference)
        aln = compress_seq['sequences']
        ref = compress_seq['reference']
        is_vcf = True
    else:
        aln = args.alignment


    # if not specified, construct default output file name with suffix _tt.nwk
    if args.output_tree:
        tree_fname = args.output_tree
    else:
        tree_fname = '.'.join(args.alignment.split('.')[:-1]) + '_tt.nwk'

    if args.timetree:
        # load meta data and covert dates to numeric
        if args.metadata is None:
            print("ERROR: meta data with dates is required for time tree reconstruction")
            return 1
        metadata, columns = read_metadata(args.metadata)
        if args.year_bounds:
            args.year_bounds.sort()
        dates = get_numerical_dates(metadata, fmt=args.date_format,
                                    min_max_year=args.year_bounds)

        # save input state string for later export
        for n in T.get_terminals():
            if n.name in metadata and 'date' in metadata[n.name]:
                n.raw_date = metadata[n.name]['date']

        if args.root and len(args.root) == 1: #if anything but a list of seqs, don't send as a list
            args.root = args.root[0]

        tt = refine(tree=T, aln=aln, ref=ref, dates=dates, confidence=args.date_confidence,
                    reroot=args.root or 'best',
                    Tc=0.01 if args.coalescent is None else args.coalescent, #use 0.01 as default coalescent time scale
                    use_marginal = args.date_inference == 'marginal',
                    branch_length_inference = args.branch_length_inference or 'auto',
                    clock_rate=args.clock_rate, clock_std=args.clock_std_dev,
                    clock_filter_iqd=args.clock_filter_iqd)

        node_data['clock'] = {'rate': tt.date2dist.clock_rate,
                              'intercept': tt.date2dist.intercept,
                              'rtt_Tmrca': -tt.date2dist.intercept/tt.date2dist.clock_rate}
        attributes.extend(['numdate', 'clock_length', 'mutation_length', 'raw_date', 'date'])
        if args.date_confidence:
            attributes.append('num_date_confidence')
    else:
        from treetime import TreeAnc
        # instantiate treetime for the sole reason to name internal nodes
        tt = TreeAnc(tree=T, aln=aln, ref=ref, gtr='JC69', verbose=1)

    node_data['nodes'] = collect_node_data(T, attributes)

    # Export refined tree and node data
    import json
    tree_success = Phylo.write(T, tree_fname, 'newick', format_branch_length='%1.8f')
    print("updated tree written to",tree_fname, file=sys.stdout)
    if args.output_node_data:
        node_data_fname = args.output_node_data
    else:
        node_data_fname = '.'.join(args.alignment.split('.')[:-1]) + '.node_data.json'

    json_success = write_json(node_data, node_data_fname)
    print("node attributes written to",node_data_fname, file=sys.stdout)

    return 0 if (tree_success and json_success) else 1
