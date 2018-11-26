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
    parser.add_argument('--pivots-per-year', type=int, default=4,
                        help="number of pivots per year")
    parser.add_argument('--minimal-frequency', type=float, default=0.05,
                        help="minimal all-time frequencies for a trajectory to be estimates")
    parser.add_argument('--output', '-o', type=str,
                        help='JSON file to save estimated frequencies to')


def run(args):

    if args.metadata is None:
        print("ERROR: meta data with dates is required for frequency estimation")
        return 1

    metadata, columns = read_metadata(args.metadata)
    dates = get_numerical_dates(metadata, fmt='%Y-%m-%d')

    dt = 1.0/args.pivots_per_year

    if args.tree:
        T = Phylo.read(args.tree, 'newick')
        # estimate tree frequencies
    elif args.alignments:
        # estimate alignment frequencies
        from .frequency_estimators import alignment_frequencies
        frequencies = {}
        for gene, fname in zip(args.gene_names, args.alignments):
            if not os.path.isfile(fname):
                print("ERROR: alignment file not found")
                return 1

            aln = MultipleSeqAlignment([seq for seq in AlignIO.read(fname, 'fasta')
                                        if not seq.name.startswith('NODE_')])
            tps = np.array([np.mean(dates[x.name]) for x in aln])
            pivots = np.arange(np.floor(tps.min()/dt)*dt, np.ceil(tps.max()/dt)*dt, dt)
            freqs = alignment_frequencies(aln, tps, pivots)
            freqs.mutation_frequencies(min_freq = args.minimal_frequency)
            frequencies.update({"%s:%d%s"%(gene, pos+1, state):list(val)
                                for (pos, state), val in freqs.frequencies.items()})

        json_success = write_json(frequencies, args.output)
        print("mutation frequencies written to", args.output, file=sys.stdout)
