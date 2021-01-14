"""
infer frequencies of mutations or clades
"""
import json, os, sys
import datetime
import numpy as np
from collections import defaultdict
from Bio import Phylo, AlignIO
from Bio.Align import MultipleSeqAlignment
import treetime.utils

from .frequency_estimators import get_pivots, alignment_frequencies, tree_frequencies
from .frequency_estimators import AlignmentKdeFrequencies, TreeKdeFrequencies, TreeKdeFrequenciesError
from .utils import read_metadata, read_node_data, write_json, get_numerical_dates


def register_arguments(parser):
    # Shared arguments
    parser.add_argument('--method', choices=["diffusion", "kde"], required=True,
                        help="method by which frequencies should be estimated")
    parser.add_argument('--metadata', type=str, required=True, metavar="FILE",
                        help="metadata including dates for given samples, as CSV or TSV")
    parser.add_argument('--regions', type=str, nargs='+', default=['global'],
                        help="region to subsample to")
    parser.add_argument("--pivot-interval", type=int, default=3,
                        help="number of units between pivots")
    parser.add_argument("--pivot-interval-units", type=str, default="months", choices=['months', 'weeks'],
                        help="space pivots by months (default) or by weeks")
    parser.add_argument('--min-date', type=numeric_date,
                        help="date to begin frequencies calculations; may be specified as an Augur-style numeric date (with the year as the integer part) or YYYY-MM-DD")
    parser.add_argument('--max-date', type=numeric_date,
                        help="date to end frequencies calculations; may be specified as an Augur-style numeric date (with the year as the integer part) or YYYY-MM-DD")

    # Tree-specific arguments
    parser.add_argument('--tree', '-t', type=str,
                        help="tree to estimate clade frequencies for")
    parser.add_argument("--include-internal-nodes", action="store_true",
                        help="calculate frequencies for internal nodes as well as tips")

    # Alignment-specific arguments
    parser.add_argument('--alignments', type=str, nargs='+',
                        help="alignments to estimate mutations frequencies for")
    parser.add_argument('--gene-names', nargs='+', type=str,
                        help="names of the sequences in the alignment, same order assumed")
    parser.add_argument('--ignore-char', type=str, default='',
                        help="character to be ignored in frequency calculations")
    parser.add_argument('--minimal-frequency', type=float, default=0.05,
                        help="minimal all-time frequencies for a trajectory to be estimates")

    # KDE-specific arguments
    parser.add_argument("--narrow-bandwidth", type=float, default=1 / 12.0, help="the bandwidth for the narrow KDE")
    parser.add_argument("--wide-bandwidth", type=float, default=3 / 12.0, help="the bandwidth for the wide KDE")
    parser.add_argument("--proportion-wide", type=float, default=0.2, help="the proportion of the wide bandwidth to use in the KDE mixture model")
    parser.add_argument("--weights", help="a dictionary of key/value mappings in JSON format used to weight KDE tip frequencies")
    parser.add_argument("--weights-attribute", help="name of the attribute on each tip whose values map to the given weights dictionary")
    parser.add_argument("--censored", action="store_true", help="calculate censored frequencies at each pivot")

    # Diffusion frequency specific arguments
    parser.add_argument('--minimal-clade-size', type=int, default=0,
                        help="minimal number of tips a clade must have for its diffusion frequencies to be reported")
    parser.add_argument('--minimal-clade-size-to-estimate', type=int, default=10,
                        help="""minimal number of tips a clade must have for its diffusion frequencies to be estimated
                                by the diffusion likelihood; all smaller clades will inherit frequencies from their
                                parents""")
    parser.add_argument("--stiffness", type=float, default=10.0, help="parameter penalizing curvature of the frequency trajectory")
    parser.add_argument("--inertia", type=float, default=0.0, help="determines how frequencies continue "
                        "in absense of data (inertia=0 -> go flat, inertia=1.0 -> continue current trend)")

    # Output arguments
    parser.add_argument('--output-format', default='auspice', choices=['auspice', 'nextflu'],
                        help="format to export frequencies JSON depending on the viewing interface")
    parser.add_argument('--output', '-o', type=str,
                        help='JSON file to save estimated frequencies to')


def format_frequencies(freq):
    return [round(x,6) for x in freq]


def run(args):
    metadata, columns = read_metadata(args.metadata)
    dates = get_numerical_dates(metadata, fmt='%Y-%m-%d')
    stiffness = args.stiffness
    inertia = args.inertia

    if args.method == "kde":
        # Load weights if they have been provided.
        if args.weights:
            with open(args.weights, "r", encoding='utf-8') as fh:
                weights = json.load(fh)

            weights_attribute = args.weights_attribute
        else:
            weights = None
            weights_attribute = None

    if args.tree:
        tree = Phylo.read(args.tree, 'newick')
        tps = []
        for tip in tree.get_terminals():
            tip.attr = {"num_date": np.mean(dates[tip.name])}
            tps.append(tip.attr["num_date"])

            # Annotate tips with metadata to enable filtering and weighting of
            # frequencies by metadata attributes.
            for key, value in metadata[tip.name].items():
                tip.attr[key] = value

        if args.method == "diffusion":
            # estimate tree frequencies
            pivots = get_pivots(tps, args.pivot_interval, args.min_date, args.max_date, args.pivot_interval_units)
            frequency_dict = {"pivots":format_frequencies(pivots)}
            frequency_dict["counts"] = {}

            for region in args.regions:
                # Omit strains sampled prior to the first pivot from frequency calculations.
                # (these tend to be reference strains included for phylogenetic context)
                if region=='global':
                    node_filter_func = lambda node: node.attr["num_date"] >= pivots[0]
                else:
                    node_filter_func = lambda node: (node.attr["region"] == region
                                                    and node.attr["num_date"] >= pivots[0])

                tree_freqs = tree_frequencies(tree, pivots, method='SLSQP',
                                              node_filter = node_filter_func,
                                              ws = max(2, tree.count_terminals()//10),
                                              stiffness = stiffness, inertia=inertia,
                                              min_clades=args.minimal_clade_size_to_estimate)

                tree_freqs.estimate_clade_frequencies()

                frequency_dict["counts"][region] = [int(x) for x in tree_freqs.counts]
                if args.output_format == "nextflu":
                    # Export frequencies in nextflu-format by region and clade id.
                    for clade_id, clade_frequencies in tree_freqs.frequencies.items():
                        frequency_dict["%s_clade:%d" % (region, clade_id)] = format_frequencies(clade_frequencies)

                else:
                    # Export frequencies in auspice-format by strain name.
                    for node in tree.find_clades(order='postorder'):
                        if node.is_terminal():
                            node.tipcount=1
                        else:
                            node.tipcount = np.sum([c.tipcount for c in node])

                        if (node.is_terminal() or args.include_internal_nodes) and node.tipcount>args.minimal_clade_size:
                            if node.name not in frequency_dict:
                                frequency_dict[node.name] = {}
                            frequency_dict[node.name][region] = format_frequencies(tree_freqs.frequencies[node.clade])

        elif args.method == "kde":
            if args.output_format == "nextflu":
                print("ERROR: nextflu format is not supported for KDE frequencies", file=sys.stderr)
                return 1

            # Estimate frequencies.
            kde_frequencies = TreeKdeFrequencies(
                sigma_narrow=args.narrow_bandwidth,
                sigma_wide=args.wide_bandwidth,
                proportion_wide=args.proportion_wide,
                pivot_frequency=args.pivot_interval,
                start_date=args.min_date,
                end_date=args.max_date,
                pivot_interval_units=args.pivot_interval_units,
                weights=weights,
                weights_attribute=weights_attribute,
                include_internal_nodes=args.include_internal_nodes,
                censored=args.censored
            )

            try:
                frequencies = kde_frequencies.estimate(tree)
            except TreeKdeFrequenciesError as e:
                print("ERROR: %s" % str(e), file=sys.stderr)
                return 1

            # Export frequencies in auspice-format by strain name.
            frequency_dict = {"pivots": list(kde_frequencies.pivots)}
            for node_name in frequencies:
                frequency_dict[node_name] = {
                    "frequencies": format_frequencies(frequencies[node_name])
                }

        write_json(frequency_dict, args.output)
        print("tree frequencies written to", args.output, file=sys.stdout)
    elif args.alignments:
        frequencies = None
        for gene, fname in zip(args.gene_names, args.alignments):
            if not os.path.isfile(fname):
                print("ERROR: alignment file not found", file=sys.stderr)
                return 1

            aln = MultipleSeqAlignment([seq for seq in AlignIO.read(fname, 'fasta')
                                        if not seq.name.startswith('NODE_')])
            tps = np.array([np.mean(dates[seq.name]) for seq in aln])

            if frequencies is None:
                pivots = get_pivots(tps, args.pivot_interval, args.min_date, args.max_date, args.pivot_interval_units)
                frequencies = {"pivots":format_frequencies(pivots)}

            if args.method == "kde":
                kde_frequencies = AlignmentKdeFrequencies(
                    sigma_narrow=args.narrow_bandwidth,
                    sigma_wide=args.wide_bandwidth,
                    proportion_wide=args.proportion_wide,
                    pivot_frequency=args.pivot_interval,
                    start_date=args.min_date,
                    end_date=args.max_date,
                    pivot_interval_units=args.pivot_interval_units,
                    weights=weights,
                    weights_attribute=weights_attribute,
                    include_internal_nodes=args.include_internal_nodes,
                    censored=args.censored
                )
                kde_frequencies.estimate(
                    aln,
                    tps
                )

                for mutation, mutation_frequencies in kde_frequencies.frequencies.items():
                    position, state = mutation.split(":")
                    frequencies["%s:%s%s" % (gene, position, state)] = format_frequencies(mutation_frequencies)
            else:
                freqs = alignment_frequencies(aln, tps, pivots, stiffness=stiffness, inertia=inertia, method='SLSQP', dtps=2.0)
                freqs.mutation_frequencies(min_freq = args.minimal_frequency, ignore_char=args.ignore_char)
                frequencies.update({"%s:%d%s" % (gene, pos+1, state): format_frequencies(mutation_frequencies)
                                    for (pos, state), mutation_frequencies in freqs.frequencies.items()})
                frequencies["%s:counts" % gene] = [int(observations_per_pivot)
                                                   for observations_per_pivot in freqs.counts]

        write_json(frequencies, args.output)
        print("mutation frequencies written to", args.output, file=sys.stdout)


def numeric_date(date):
    """
    Converts the given *date* string to a :py:class:`float`.
    *date* may be given as a number (a float) with year as the integer part, or
    in the YYYY-MM-DD (ISO 8601) syntax.
    >>> numeric_date("2020.42")
    2020.42
    >>> numeric_date("2020-06-04")
    2020.42486...
    """
    try:
        return float(date)
    except ValueError:
        return treetime.utils.numeric_date(datetime.date(*map(int, date.split("-", 2))))
