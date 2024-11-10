"""
infer frequencies of mutations or clades
"""
import json, os, sys
import numpy as np
from Bio import Phylo, AlignIO
from Bio.Align import MultipleSeqAlignment

from .argparse_ import ExtendOverwriteDefault
from .errors import AugurError
from .frequency_estimators import get_pivots, alignment_frequencies, tree_frequencies
from .frequency_estimators import AlignmentKdeFrequencies, TreeKdeFrequencies, TreeKdeFrequenciesError
from .dates import numeric_date_type, SUPPORTED_DATE_HELP_TEXT, get_numerical_dates
from .io.file import open_file
from .io.metadata import DEFAULT_DELIMITERS, DEFAULT_ID_COLUMNS, METADATA_DATE_COLUMN, InvalidDelimiter, Metadata, read_metadata
from .utils import write_json

REGION_COLUMN = 'region'
DEFAULT_REGION = 'global'

def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("frequencies", help=__doc__)
    # Shared arguments
    parser.add_argument('--method', choices=["diffusion", "kde"], required=True,
                        help="method by which frequencies should be estimated")
    parser.add_argument('--metadata', type=str, required=True, metavar="FILE",
                        help="metadata including dates for given samples")
    parser.add_argument('--metadata-delimiters', default=DEFAULT_DELIMITERS, nargs="+", action=ExtendOverwriteDefault,
                        help="delimiters to accept when reading a metadata file. Only one delimiter will be inferred.")
    parser.add_argument('--metadata-id-columns', default=DEFAULT_ID_COLUMNS, nargs="+", action=ExtendOverwriteDefault,
                        help="names of possible metadata columns containing identifier information, ordered by priority. Only one ID column will be inferred.")
    parser.add_argument('--regions', type=str, nargs='+', action=ExtendOverwriteDefault, default=[DEFAULT_REGION],
                        help="region to filter to. " \
                            f"Regions should match values in the {REGION_COLUMN!r} column of the metadata file " \
                            f"if specifying values other than the default {DEFAULT_REGION!r} region.")
    parser.add_argument("--pivot-interval", type=int, default=3,
                        help="number of units between pivots")
    parser.add_argument("--pivot-interval-units", type=str, default="months", choices=['months', 'weeks'],
                        help="space pivots by months (default) or by weeks")
    parser.add_argument('--min-date', type=numeric_date_type,
                        help=f"date to begin frequencies calculations; may be specified as: {SUPPORTED_DATE_HELP_TEXT}")
    parser.add_argument('--max-date', type=numeric_date_type,
                        help=f"date to end frequencies calculations; may be specified as: {SUPPORTED_DATE_HELP_TEXT}")

    # Tree-specific arguments
    parser.add_argument('--tree', '-t', type=str,
                        help="tree to estimate clade frequencies for")
    parser.add_argument("--include-internal-nodes", action="store_true",
                        help="calculate frequencies for internal nodes as well as tips")

    # Alignment-specific arguments
    parser.add_argument('--alignments', type=str, nargs='+', action='extend',
                        help="alignments to estimate mutations frequencies for")
    parser.add_argument('--gene-names', nargs='+', action='extend', type=str,
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
    return parser


def format_frequencies(freq):
    return [round(x,6) for x in freq]


def run(args):
    try:
        metadata_object = Metadata(args.metadata, args.metadata_delimiters, args.metadata_id_columns)
    except InvalidDelimiter:
        raise AugurError(
            f"Could not determine the delimiter of {args.metadata!r}. "
            f"Valid delimiters are: {args.metadata_delimiters!r}. "
            "This can be changed with --metadata-delimiters."
        )

    columns_to_load = [metadata_object.id_column, METADATA_DATE_COLUMN]
    if args.weights_attribute:
        columns_to_load.append(args.weights_attribute)

    filter_to_region = any(region != DEFAULT_REGION for region in args.regions)
    if filter_to_region:
        columns_to_load.append(REGION_COLUMN)

    metadata = read_metadata(
        args.metadata,
        delimiters=[metadata_object.delimiter],
        columns=columns_to_load,
        id_columns=[metadata_object.id_column],
        dtype="string",
    )
    dates = get_numerical_dates(metadata, fmt='%Y-%m-%d')
    stiffness = args.stiffness
    inertia = args.inertia

    # For KDE
    weights = None
    weights_attribute = None

    if args.method == "kde" and args.weights:
        # Load weights if they have been provided.
        with open_file(args.weights, "r") as fh:
            weights = json.load(fh)

        weights_attribute = args.weights_attribute

    if args.tree:
        tree = Phylo.read(args.tree, 'newick')
        tps = []
        for tip in tree.get_terminals():
            tip.attr = {"num_date": np.mean(dates[tip.name])}
            tps.append(tip.attr["num_date"])

            if weights_attribute:
                # Annotate tip with weight attribute.
                tip.attr[weights_attribute] = metadata.loc[tip.name, weights_attribute]

            if filter_to_region:
                tip.attr[REGION_COLUMN] = metadata.loc[tip.name, REGION_COLUMN]

        if args.method == "diffusion":
            # estimate tree frequencies
            pivots = get_pivots(tps, args.pivot_interval, args.min_date, args.max_date, args.pivot_interval_units)
            frequency_dict = {"pivots":format_frequencies(pivots)}
            frequency_dict["counts"] = {}

            for region in args.regions:
                # Omit strains sampled prior to the first pivot from frequency calculations.
                # (these tend to be reference strains included for phylogenetic context)
                if region==DEFAULT_REGION:
                    node_filter_func = lambda node: node.attr["num_date"] >= pivots[0]
                else:
                    node_filter_func = lambda node: (node.attr[REGION_COLUMN] == region
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
