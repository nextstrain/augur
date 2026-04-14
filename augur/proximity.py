"""
Find proximal sequences for focal sequences vs contextual sequences.
"""

import argparse
import math
import sys
import time
import threading
from textwrap import dedent
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
from numpy.typing import NDArray
from augur.errors import AugurError
from augur.io.print import print_err
from augur.io.sequences import read_sequences, subset_fasta
from augur.utils import nthreads_value

"""
Currently this uses hamming distances as proximity, with different options around how we treat missing (non-ATGC) data.

This builds off our ncov implementation <https://github.com/nextstrain/ncov/blob/master/scripts/get_distance_to_focal_set.py>
with some important differences. Here we compute prioirities _per focal sequence_ so that we pick contextual sequences
for every focal strain, even if there are no close contextual sequences for that strain. Secondly, our approach here is reference-
free (although inputs must be aligned).

When there's tiebreaks involved in picking k-sequences, we pick the sequences with the fewest Ns. Tiebreaks
in the number of Ns are resolved alphabetically (and thus the result is deterministic).

Future directions
=================

* Add temporal metadata so we can have criteria such as "max (nuc) distance X and max temporal distance Y".

* Hamming distance in AA-space, related to <https://github.com/nextstrain/augur/issues/820>

* Different proximity methods, such as tree-based <https://github.com/flu-crew/parnas>, or tree-informed
  (e.g. use hamming distance to pick a cloud of sequences, build a small tree, pick the closest on the tree)

* More nuance around non-ATGC characters (especially gaps)


Performance optimisations
=========================

* We currently load all (context) sequences into memory. We could instead load these one batch at a time,
  calculate neighbours for all focal seqs x this batch, and then combine all the neighbours at the end.
  At the moment, 350k context sequences of 2.5kb (influenza-like) are <2Gb memory, so I don't
  think it's worth the code complexity at the moment.

"""


N_VALUE: int = ord('n')
VALID_NUCS: set[int] = {ord('a'), ord('t'), ord('c'), ord('g')} # anything else will be considered "n"

proximity_argument_descriptions: dict[str, str] = {
    'method': "Proximity approach used",
    'context_sequences': "FASTA file with aligned contextual sequences",
    'focal_sequences': "FASTA file with aligned focal sequences to find neighbors for",
    'output_strains': "output file with one neighbor strain name per line",
    'output_matches': "optional TSV file with columns: focal strain, context strain, distance",
    'output_sequences': "All proximal strains found",
    'k': "number of nearest neighbors to find per focal strain",
    'max_distance': "maximum distance threshold for considering a sequence to match",
    'no_progress': "Don't print ongoing progress output",
    'ignore_missing_data': dedent("""\
        All non-ATGC bases are converted to 'N', and then:
        - 'none' treats 'N' as a normal base for comparison purposes;
        - 'all' ignores positions where either sequence is N;
        - 'flanking' ignores runs of Ns at the start/end of each sequence."""),
    'nthreads': "Number of threads to use for parallel processing. Use 'auto' to use all available cores.",
}
"""
`augur proximity` argument descriptions, stored as a dict so we can re-use in `augur subsample` related code
"""

def register_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument('--method', default="hamming", choices=['hamming'], help=proximity_argument_descriptions['method'])
    parser.add_argument('--context-sequences', required=True, help=proximity_argument_descriptions['context_sequences'])
    parser.add_argument('--focal-sequences', required=True, help=proximity_argument_descriptions['focal_sequences'])
    parser.add_argument('--output-strains', required=True, help=proximity_argument_descriptions['output_strains'])
    parser.add_argument('--output-matches', help=proximity_argument_descriptions['output_matches'])
    parser.add_argument('--output-sequences', metavar="FASTA", help=proximity_argument_descriptions['output_sequences'])
    parser.add_argument('--k', type=int, default=5, help=proximity_argument_descriptions['k'])
    parser.add_argument('--max-distance', type=int, default=4, help=proximity_argument_descriptions['max_distance'])
    parser.add_argument('--no-progress', action="store_true", default=False, help=proximity_argument_descriptions['no_progress'])
    parser.add_argument('--ignore-missing-data', choices=['none', 'all', 'flanking'], default='none', help=proximity_argument_descriptions['ignore_missing_data'])
    parser.add_argument('--nthreads', type=nthreads_value, default=1, help=proximity_argument_descriptions['nthreads'])


def register_parser(parent_subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = parent_subparsers.add_parser("proximity", help=__doc__)
    register_arguments(parser)
    return parser


def to_numpy_array(seq: str) -> NDArray[np.int8]:
    """Convert a nucleotide sequence string to a numpy integer array.

    Uses raw byte values of the lowercase string (a=97, c=99, g=103, t=116).
    Non-ATCG characters are replaced with the N sentinel value (n=110).
    """
    arr = np.frombuffer(str(seq).lower().encode('utf-8'), dtype=np.int8).copy()
    mask = np.ones(len(arr), dtype=bool)
    for nuc in VALID_NUCS:
        mask &= (arr != nuc)
    arr[mask] = N_VALUE
    return arr


def get_valid_range(seq: NDArray[np.int8]) -> tuple[int, int]:
    """Return (start, end) indices of the non-flanking-N region of *seq*.
    start is the (0-based) index of the first non-N character
    end is the (0-based) index of the first N character of the trailing run
    """
    n = len(seq)
    start = 0
    while start < n and seq[start] == N_VALUE:
        start += 1
    end = n
    while end > start and seq[end - 1] == N_VALUE:
        end -= 1
    return (start, end)


def select_top_k(
        all_distances: NDArray[np.int32], # shape: (n_context, 1)
        context_names: list[str],
        context_matrix: NDArray[np.int8],
        k: int,
        max_distance: int) -> list[dict[str, str | int]]:
    """
    Select top-k sequences (within the max_distance threshold) from a vector of
    distance scores. Results are sorted by distance counts, then alphabetically
    for those with the same distance.
    
    To resolve tiebreaks (e.g. for k=3, max_distance=5, distances=[1,1,2,2,2,...],
    which of the distance=2 strains do we take?) we choose the strains with the
    fewest "N"s, with remaining ties resolved alpabetically.
    """
    # get the (context) indexes which are candidates for selection
    within_threshold = np.where(all_distances <= max_distance)[0] # array
    if len(within_threshold) == 0:
        return []
    if len(within_threshold) <= k:
        selected = sorted(within_threshold, key=lambda i: (all_distances[i], context_names[i]))
        return [{"strain": context_names[i], "distance": int(all_distances[i])} for i in selected]

    # More than k candidates - work out which distance value straddles the k-boundary
    # I.e. distance<cutoff we'd get fewer than k, distance<=cutoff we'd get more than k.
    cutoff = int(np.partition(all_distances[within_threshold], k - 1)[k - 1])
    threshold_distances = all_distances[within_threshold]
    below = within_threshold[threshold_distances < cutoff]
    at_cutoff = within_threshold[threshold_distances == cutoff]
    need = k - len(below) # num strains at the cutoff boundary we need to choose
    # Only compute counts of "N" for sequences at the cutoff distance (for performance)
    at_cutoff_n_counts = np.sum(context_matrix[at_cutoff] == N_VALUE, axis=1)
    at_cutoff_sorted = sorted(
        range(len(at_cutoff)),
        key=lambda j: (int(at_cutoff_n_counts[j]), context_names[at_cutoff[j]]),
    )
    at_cutoff_selected = [at_cutoff[at_cutoff_sorted[j]] for j in range(need)]

    below_sorted = sorted(below, key=lambda i: (all_distances[i], context_names[i]))
    selected = list(below_sorted) + at_cutoff_selected
    return [{"strain": context_names[i], "distance": int(all_distances[i])} for i in selected]

def distance_fn(
        ignore_missing_data: str,
        focal_seq: NDArray[np.int8],
        context_matrix: NDArray[np.int8],
        focal_valid_range: tuple[int, int] | None = None,
        context_valid_mask: NDArray[np.bool_] | None = None,
        ) -> Callable[[NDArray[np.int8], int, int], NDArray[np.int8]]:
    """
    Creates a function to be applied to a batch of context sequences
    which computes distances of this focal_seq to the batch of context sequences.
    
    Returned function: (batch, batch_start, batch_end) -> distances
    
    How we count "N"s depends on the *ignore_missing_data*
    """
    if ignore_missing_data=='none':
        # "N" isn't special, just count all differences
        def _compute_batch_none(batch, _batch_start, _batch_end):
            return np.sum(focal_seq != batch, axis=1)
        return _compute_batch_none
    if ignore_missing_data=='flanking':
        assert focal_valid_range is not None and context_valid_mask is not None
        q_start, q_end = focal_valid_range
        focal_mask = np.zeros(context_matrix.shape[1], dtype=bool) # all values `False`
        focal_mask[q_start:q_end] = True
        # Compute the combined mask per-batch (not upfront for all contexts) to avoid
        # allocating a full n_context × seq_len array per focal seq, which causes cache
        # thrashing and poor thread scaling.
        def _compute_batch_flanking(batch, batch_start, batch_end):
            batch_mask = context_valid_mask[batch_start:batch_end] & focal_mask
            return np.sum((focal_seq != batch) & batch_mask, axis=1)
        return _compute_batch_flanking
    if ignore_missing_data=='all':
        focal_valid = focal_seq != N_VALUE
        def _compute_batch_all(batch, _batch_start, _batch_end):
            valid = focal_valid & (batch != N_VALUE) # True if both context & focal seqs are non-"N"
            return np.sum((focal_seq != batch) & valid, axis=1)
        return _compute_batch_all
    raise AugurError(f"Unexpcted {ignore_missing_data:=}")


def load_focal_sequences(fname: str) -> tuple[dict[str, NDArray[np.int8]], int]:
    """
    Load all focal sequences into memory as numpy arrays
    """
    # We could batch these pretty easily for memory reasons, but since
    # focal << context we qould be better to focus on batching context for memory
    seq_len:int = 0
    focal: dict[str, NDArray[np.int8]] = {}
    for record in read_sequences(fname):
        q = to_numpy_array(str(record.seq))
        if not seq_len:
            seq_len = len(q)
        elif seq_len!=len(q):
            raise AugurError(dedent(f"""\
                When reading focal samples, {record.name!r} is {len(q):,}nt, which differs
                from previous focal sequences ({seq_len:,}nt). Focal sequences must be aligned
                for proximity calculations.
                """))
        focal[record.name] = q
    if not focal:
        raise AugurError(f"No focal sequences found in '{fname}'")
    return focal, seq_len

def load_context(fname: str, skip_strains: set[str], seq_len: int) \
        -> tuple[NDArray[np.int8], list[str], int]:
    """
    Load context sequences directly into a 2D numpy matrix
    (avoids the need for intermediate structures which create memory bottlenecks).
    This means we need to know the number of sequences ahead of time, so we
    start with a guess and then expand the matrix as needed.

    Returns (context_matrix, context_names, skip_count).
    """
    INITIAL_CAPACITY = 32_768 # 2^15
    context_matrix = np.empty((INITIAL_CAPACITY, seq_len), dtype=np.int8)
    context_names: list[str] = []
    skip_count = 0
    n = 0
    for record in read_sequences(fname):
        if record.name in skip_strains:
            skip_count += 1
            continue
        arr = to_numpy_array(str(record.seq))
        if len(arr) != seq_len:
            raise AugurError(dedent(f"""\
                When reading sequences, {record.name!r} has length {len(arr):,}nt, which differs
                from previous focal sequences ({seq_len:,}nt). All sequences must be aligned and
                the same length for proximity calculations.
                """))
        if n >= context_matrix.shape[0]:
            context_matrix = np.resize(context_matrix, (context_matrix.shape[0] * 2, seq_len))
        context_matrix[n] = arr
        context_names.append(record.name)
        n += 1
    if n == 0:
        raise AugurError(f"No sequences found in context sequences '{fname}'")
    # trim to actual size
    # Note: this line triggered a mypy error in CI but not locally
    context_matrix = context_matrix[:n]  # type: ignore[assignment]
    return context_matrix, context_names, skip_count

def flanking_masks(
        focal: dict[str, NDArray[np.int8]],
        context_matrix: NDArray[np.int8],
        ) -> tuple[dict[str, tuple[int, int]], NDArray[np.bool_]]:
    """
    Compute ranges of flanking Ns for each focal sequence, and the same for
    contextual sequences in the form of a 2d mask matrix the same size as context_matrix
    """
    focal_valid_ranges = {name: get_valid_range(seq) for name, seq in focal.items()}
    context_valid_mask = np.ones_like(context_matrix, dtype=bool)
    for i in range(context_matrix.shape[0]):
        # Note: could vectorise if this becomes a bottleneck, but in real-world
        # influenza testing this takes <5% of the runtime.
        start, end = get_valid_range(context_matrix[i])
        context_valid_mask[i, :start] = False
        context_valid_mask[i, end:] = False
    return focal_valid_ranges, context_valid_mask

def _make_progress_callback(total: int, print_partial_progress: bool) -> Callable[[], None]:
    """Return a thread-safe callback that prints progress to stderr.
    The running count of focal seqs processed is stored here, which means it must be
    called every time a focal sequence is processed for the stats to be valid.
    We could throttle the actual `print` calls in the callback if the output's
    too verbose.
    """
    lock = threading.Lock()
    state = {'completed': 0, 'start_time': time.time()}

    def callback():
        with lock:
            state['completed'] += 1
            completed = state['completed']
            elapsed = time.time() - state['start_time']
            eta = elapsed / completed * (total - completed) if completed > 0 else 0
            pct = completed / total * 100
            if print_partial_progress:
                print(f"\rProgress: {completed}/{total} ({pct:.0f}%), "
                    f"Elapsed: {elapsed:.1f}s, ETA: {eta:.1f}s       ",
                    end='', file=sys.stderr, flush=True)
            if completed == total:
                print(f"\nProximity calculations complete. Total time: {elapsed:.1f}s", file=sys.stderr, flush=True)

    return callback


def _process_focal(
    focal_items: list[tuple[str, NDArray[np.int8]]],
    context_matrix: NDArray[np.int8],
    context_names: list[str],
    k: int,
    max_distance: int,
    ignore_missing_data: str,
    focal_valid_ranges: dict[str, tuple[int, int]] | None,
    context_valid_mask: NDArray[np.bool_] | None,
    on_focal_done: Callable[[], None] | None = None,
    cancel_event: threading.Event | None = None,
) -> dict[str, list[dict[str, str | int]]]:
    """
    Process a chunk of focal sequences, returning {focal_strain: [{"strain": name, "distance": d}, ...]}.
    (This function is intended for parallalisation - i.e. run this function multiple times
    on different threads.)

    Computes hamming distances from focal sequences to all rows of *context_matrix*
    in batches to control peak memory usage. Returns the *k* nearest neighbors.

    *ignore_missing_data* controls N handling ("none", "all", or "flanking").
    For "flanking", *focal_valid_range* and *context_valid_mask* must be provided.
    """
    results = {}
    for focal_strain, focal_seq in focal_items:
        if cancel_event and cancel_event.is_set():
            break
        focal_valid_range = focal_valid_ranges[focal_strain] if focal_valid_ranges else None
        compute_batch = distance_fn(ignore_missing_data, focal_seq, context_matrix, focal_valid_range, context_valid_mask)

        n_context = context_matrix.shape[0]
        BATCH_SIZE = 10_000 # (tunable for memory reasons if needed)
        all_distances = np.empty(n_context, dtype=np.int32)

        for batch_start in range(0, n_context, BATCH_SIZE):
            batch_end = min(batch_start + BATCH_SIZE, n_context)
            batch = context_matrix[batch_start:batch_end]
            all_distances[batch_start:batch_end] = compute_batch(batch, batch_start, batch_end)

        results[focal_strain] = select_top_k(all_distances, context_names, context_matrix, k, max_distance)
        if on_focal_done:
            on_focal_done()
    return results


def run(args: argparse.Namespace) -> None:

    focal, seq_len = load_focal_sequences(args.focal_sequences)
    print_err(f"Read {len(focal):,} focal sequences from {args.focal_sequences}")
    
    context_matrix, context_names, sequences_skip_count = load_context(args.context_sequences, set(focal.keys()), seq_len)
    print_err(f"Loaded {len(context_names):,} comparison sequences from {args.context_sequences}. Excluded {sequences_skip_count} as they were in the focal sequences.")

    # Precompute valid ranges / mask for flanking mode
    ignore_missing_data: str = args.ignore_missing_data
    focal_valid_ranges: dict[str, tuple[int, int]] | None = None
    context_valid_mask: NDArray[np.bool_] | None = None
    if ignore_missing_data == 'flanking':
        focal_valid_ranges, context_valid_mask = flanking_masks(focal, context_matrix)
    
    # Find neighbors for each focal strain
    focal_items = list(focal.items())
    nthreads = min(args.nthreads, len(focal_items))  # don't spawn more threads than focal seqs

    cancel_event = threading.Event()

    _process_focal_kwargs = dict(
        context_matrix=context_matrix,
        context_names=context_names,
        k=args.k,
        max_distance=args.max_distance,
        ignore_missing_data=ignore_missing_data,
        focal_valid_ranges=focal_valid_ranges,
        context_valid_mask=context_valid_mask,
        on_focal_done=_make_progress_callback(len(focal_items), not args.no_progress),
        cancel_event=cancel_event,
    )

    if nthreads <= 1:
        # Fast path: no threading overhead
        results = _process_focal(focal_items, **_process_focal_kwargs)
    else:
        # Chunk focal seqs across threads
        chunk_size = math.ceil(len(focal_items) / nthreads)
        chunks = [focal_items[i:i + chunk_size] for i in range(0, len(focal_items), chunk_size)]
        print_err(f"Processing {len(focal_items):,} focal sequences across {len(chunks):,} threads")

        results = {}
        with ThreadPoolExecutor(max_workers=nthreads) as executor:
            futures = [executor.submit(_process_focal, chunk, **_process_focal_kwargs) for chunk in chunks]
            try:
                for future in as_completed(futures):
                    results.update(future.result())
            except BaseException: # BaseException captures ctrl-c
                cancel_event.set()
                executor.shutdown(wait=False, cancel_futures=True)
                raise

    all_neighbour_strain_names: set[str] = {  # pyright: ignore[reportAssignmentType]
        match["strain"]  # type: ignore[misc]
        for matches in results.values()
        for match in matches
    }
    print_err(f"Found {len(all_neighbour_strain_names):,} unique neighbor strains for {len(focal)} focal strains")

    # Write output in deterministic (sorted) order
    with open(args.output_strains, 'w') as fh:
        for strain_name in sorted(all_neighbour_strain_names):
            print(strain_name, file=fh)
    print_err(f"Wrote {len(all_neighbour_strain_names):,} strains to {args.output_strains}")

    # Write proximal sequences FASTA if requested
    if args.output_sequences:
        subset_fasta(args.context_sequences, args.output_sequences, args.output_strains, args.nthreads)
        print_err(f"Wrote proximal sequences to {args.output_sequences}")

    # Write per-focal matches TSV if requested
    if args.output_matches:
        with open(args.output_matches, 'w') as fh:
            print("focal_strain\tcontext_strain\tdistance", file=fh)
            for focal_strain, matches in results.items():
                for match in matches:
                    print(f"{focal_strain}\t{match['strain']}\t{match['distance']}", file=fh)
        print_err(f"Wrote matches to {args.output_matches}")
