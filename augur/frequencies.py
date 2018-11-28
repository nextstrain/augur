"""
infer frequencies of mutations or clades
"""
import json, os, sys
import numpy as np
from collections import defaultdict
from Bio import Phylo, AlignIO
from Bio.Align import MultipleSeqAlignment
from .utils import read_metadata, read_node_data, write_json, get_numerical_dates


def register_arguments(parser):
    parser.add_argument('--tree', '-t', type=str,
                        help="tree to estimate clade frequencies for")
    parser.add_argument('--alignments', type=str, nargs='+',
                        help="alignments to estimate mutations frequencies for")
    parser.add_argument('--metadata', type=str,
                        help="alignment to estimate mutations frequencies for")
    parser.add_argument('--gene-names', nargs='+', type=str,
                        help="names of the sequences in the alignment, same order assumed")
    parser.add_argument('--regions', type=str, nargs='+', default=['global'],
                        help="region to subsample to")
    parser.add_argument('--min-date', type=float,
                        help="minimal pivot value")
    parser.add_argument('--pivots-per-year', type=int, default=4,
                        help="number of pivots per year")
    parser.add_argument('--minimal-frequency', type=float, default=0.05,
                        help="minimal all-time frequencies for a trajectory to be estimates")
    parser.add_argument('--output', '-o', type=str,
                        help='JSON file to save estimated frequencies to')


def format_frequencies(freq):
    return list(freq)

def run(args):

    if args.metadata is None:
        print("ERROR: meta data with dates is required for frequency estimation")
        return 1

    metadata, columns = read_metadata(args.metadata)
    dates = get_numerical_dates(metadata, fmt='%Y-%m-%d')
    dt = 1.0/args.pivots_per_year
    stiffness = 5.0

    if args.tree:
        from .frequency_estimators import tree_frequencies
        tree = Phylo.read(args.tree, 'newick')
        tps = []
        for x in tree.get_terminals():
            x.num_date = np.mean(dates[x.name])
            tps.append(x.num_date)
            x.region = metadata[x.name]["region"]

        first_pivot = args.min_date if args.min_date else np.floor(np.min(tps)/dt)*dt
        pivots = np.arange(first_pivot, np.ceil(np.max(tps)/dt)*dt, dt)
        frequency_dict = {"pivots":format_frequencies(pivots)}

        # estimate tree frequencies
        # Omit strains sampled prior to the first pivot from frequency calculations.
        for region in args.regions:
            if region=='global':
                node_filter_func = lambda node: node.num_date >= first_pivot
            else:
                node_filter_func = lambda node: (node.region == region
                                                and node.num_date >= first_pivot)

            tree_freqs = tree_frequencies(tree, pivots, method='SLSQP',
                                          node_filter = node_filter_func,
                                          ws = max(2, tree.count_terminals()//10),
                                          stiffness = stiffness)

            tree_freqs.estimate_clade_frequencies()
            for x,val in tree_freqs.frequencies.items():
                frequency_dict["%s_clade:%d"%(region,x)] = format_frequencies(val)

        json_success = write_json(frequency_dict, args.output)
        print("tree frequencies written to", args.output, file=sys.stdout)

    elif args.alignments:
        from .frequency_estimators import alignment_frequencies
        frequencies = None
        for gene, fname in zip(args.gene_names, args.alignments):
            if not os.path.isfile(fname):
                print("ERROR: alignment file not found")
                return 1

            aln = MultipleSeqAlignment([seq for seq in AlignIO.read(fname, 'fasta')
                                        if not seq.name.startswith('NODE_')])
            tps = np.array([np.mean(dates[x.name]) for x in aln])
            if frequencies is None:
                pivots = np.arange(np.floor(tps.min()/dt)*dt, np.ceil(tps.max()/dt)*dt, dt)
                frequencies = {"pivots":format_frequencies(pivots)}

            freqs = alignment_frequencies(aln, tps, pivots)
            freqs.mutation_frequencies(min_freq = args.minimal_frequency)
            frequencies.update({"%s:%d%s"%(gene, pos+1, state):format_frequencies(val)
                                for (pos, state), val in freqs.frequencies.items()})

        json_success = write_json(frequencies, args.output)
        print("mutation frequencies written to", args.output, file=sys.stdout)
