"""
Find proximal sequences for query strains vs contextual sequences.
"""

import argparse
import math
import sys
import time
import threading
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
from numpy.typing import NDArray
from augur.errors import AugurError
from augur.io.print import print_err
from augur.io.sequences import read_sequences
from augur.utils import nthreads_value

"""
Currently this uses hamming distances as proximity, with different options around how we treat missing (non-ATGC) data.

This builds off our ncov implementation <https://github.com/nextstrain/ncov/blob/master/scripts/get_distance_to_focal_set.py>
with some important differences. Here we compute prioirities _per query (focal) sequence_ so that we pick contextual sequences
for every query strain, even if there are no close contextual sequences for that strain. Secondly, our approach here is reference-
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
  calculate neighbours for all queries x this batch, and then combine all the neighbours at the end.
  At the moment, 350k context sequences of 2.5kb (influenza-like) are <2Gb memory, so I don't
  think it's worth the code complexity at the moment.

"""


N_VALUE: int = ord('n')
VALID_NUCS: set[int] = {ord('a'), ord('t'), ord('c'), ord('g')} # anything else will be considered "n"


def register_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument('--method', default="hamming", choices=['hamming'], help="Proximity approach used")
    parser.add_argument('--sequences', required=True, help="FASTA file with aligned contextual sequences")
    parser.add_argument('--query', required=True, help="FASTA file with aligned query sequences to find neighbors for")
    parser.add_argument('--output-strains', required=True, help="output file with one neighbor strain name per line")
    parser.add_argument('--output-matches', help="optional TSV file with columns: query strain, context strain, distance")
    parser.add_argument('--k', type=int, default=5, help="number of nearest neighbors to find per query strain")
    parser.add_argument('--max-distance', type=int, default=4, help="maximum distance threshold for considering a sequence to match")
    parser.add_argument('--no-progress', action="store_true", default=False, help="Dont print ongoing progress output")
    parser.add_argument('--missing-data', choices=['none', 'all', 'flanking'], default='none',
                        help="All non-ATGC bases are converted to 'N', and then: "
                             "'none' treats 'N' as a normal base for comparison purposes; "
                             "'all' ignores positions where either sequence is N; "
                             "'flanking' ignores runs of Ns at the start/end of each sequence.")
    parser.add_argument('--nthreads', type=nthreads_value, default=1,
                        help="Number of threads to use for parallel query processing. "
                             "Use 'auto' to use all available cores.")


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
        missing_data: str,
        query_seq: NDArray[np.int8],
        context_matrix: NDArray[np.int8],
        query_valid_range: tuple[int, int] | None = None,
        context_valid_mask: NDArray[np.bool_] | None = None,
        ) -> Callable[[NDArray[np.int8], int, int], NDArray[np.int8]]:
    """
    Creates a function to be applied to a batch of context sequences
    which computes distances of this query_seq to the batch of context sequences.
    
    Returned function: (batch, batch_start, batch_end) -> distances
    
    How we count "N"s depends on *missing_data*
    """
    if missing_data=='none':
        # "N" isn't special, just count all differences
        def _compute_batch_none(batch, _batch_start, _batch_end):
            return np.sum(query_seq != batch, axis=1)
        return _compute_batch_none
    if missing_data=='flanking':
        assert query_valid_range is not None and context_valid_mask is not None
        q_start, q_end = query_valid_range
        query_mask = np.zeros(context_matrix.shape[1], dtype=bool) # all values `False`
        query_mask[q_start:q_end] = True
        # Compute the combined mask per-batch (not upfront for all contexts) to avoid
        # allocating a full n_context × seq_len array per query, which causes cache
        # thrashing and poor thread scaling.
        def _compute_batch_flanking(batch, batch_start, batch_end):
            batch_mask = context_valid_mask[batch_start:batch_end] & query_mask
            return np.sum((query_seq != batch) & batch_mask, axis=1)
        return _compute_batch_flanking
    if missing_data=='all':
        query_valid = query_seq != N_VALUE
        def _compute_batch_all(batch, _batch_start, _batch_end):
            valid = query_valid & (batch != N_VALUE) # True if both context & query are non-"N"
            return np.sum((query_seq != batch) & valid, axis=1)
        return _compute_batch_all
    raise AugurError(f"Unexpcted {missing_data:=}")


def load_query(fname: str) -> tuple[dict[str, NDArray[np.int8]], int]:
    """
    Load all query sequences into memory as numpy arrays
    """
    # We could batch these pretty easily for memory reasons, but since
    # queries << context we qould be better to focus on batching context for memory
    seq_len:int = 0
    queries: dict[str, NDArray[np.int8]] = {}
    for record in read_sequences(fname):
        q = to_numpy_array(str(record.seq))
        if not seq_len:
            seq_len = len(q)
        elif seq_len!=len(q):
            raise AugurError(f"Query sequence '{record.name}' has length {len(q)} "
                f"but previous queries had length {seq_len}")
        queries[record.name] = q
    if not queries:
        raise AugurError(f"No query sequences found in '{fname}'")
    return queries, seq_len

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
            raise AugurError(f"Context sequence '{record.name}' has length {len(arr)} "
                             f"but expected length {seq_len}")
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
        queries: dict[str, NDArray[np.int8]],
        context_matrix: NDArray[np.int8],
        ) -> tuple[dict[str, tuple[int, int]], NDArray[np.bool_]]:
    """
    Compute ranges of flanking Ns for each query sequence, and the same for
    contextual sequences in the form of a 2d mask matrix the same size as context_matrix
    """
    query_valid_ranges = {name: get_valid_range(seq) for name, seq in queries.items()}
    context_valid_mask = np.ones_like(context_matrix, dtype=bool)
    for i in range(context_matrix.shape[0]):
        # Note: could vectorise if this becomes a bottleneck, but in real-world
        # influenza testing this takes <5% of the runtime.
        start, end = get_valid_range(context_matrix[i])
        context_valid_mask[i, :start] = False
        context_valid_mask[i, end:] = False
    return query_valid_ranges, context_valid_mask

def _make_progress_callback(total: int, print_partial_progress: bool) -> Callable[[], None]:
    """Return a thread-safe callback that prints progress to stderr.
    The running count of queries processed is stored here, which means it must be
    called every time a query is processed for the stats to be valid.
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


def _process_queries(
    query_items: list[tuple[str, NDArray[np.int8]]],
    context_matrix: NDArray[np.int8],
    context_names: list[str],
    k: int,
    max_distance: int,
    missing_data: str,
    query_valid_ranges: dict[str, tuple[int, int]] | None,
    context_valid_mask: NDArray[np.bool_] | None,
    on_query_done: Callable[[], None] | None = None,
    cancel_event: threading.Event | None = None,
) -> dict[str, list[dict[str, str | int]]]:
    """
    Process a chunk of query sequences, returning {query_strain: [{"strain": name, "distance": d}, ...]}.
    (This function is intended for parallalisation - i.e. run this function multiple times
    on different threads.)

    Computes hamming distances from *query_seq* to all rows of *context_matrix*
    in batches to control peak memory usage. Returns the *k* nearest neighbors.

    *missing_data* controls N handling ("none", "all", or "flanking").
    For "flanking", *query_valid_range* and *context_valid_mask* must be provided.
    """
    results = {}
    for query_strain, query_seq in query_items:
        if cancel_event and cancel_event.is_set():
            break
        query_valid_range = query_valid_ranges[query_strain] if query_valid_ranges else None
        compute_batch = distance_fn(missing_data, query_seq, context_matrix, query_valid_range, context_valid_mask)

        n_context = context_matrix.shape[0]
        BATCH_SIZE = 10_000 # (tunable for memory reasons if needed)
        all_distances = np.empty(n_context, dtype=np.int32)

        for batch_start in range(0, n_context, BATCH_SIZE):
            batch_end = min(batch_start + BATCH_SIZE, n_context)
            batch = context_matrix[batch_start:batch_end]
            all_distances[batch_start:batch_end] = compute_batch(batch, batch_start, batch_end)

        results[query_strain] = select_top_k(all_distances, context_names, context_matrix, k, max_distance)
        if on_query_done:
            on_query_done()
    return results


def run(args: argparse.Namespace) -> None:

    queries, seq_len = load_query(args.query)
    print_err(f"Read {len(queries):,} query sequences from {args.query}")
    
    context_matrix, context_names, sequences_skip_count = load_context(args.sequences, set(queries.keys()), seq_len)
    print_err(f"Loaded {len(context_names):,} comparison sequences from {args.sequences}. Excluded {sequences_skip_count} as they were in the focal sequences.")

    # Precompute valid ranges / mask for flanking mode
    missing_data: str = args.missing_data
    query_valid_ranges: dict[str, tuple[int, int]] | None = None
    context_valid_mask: NDArray[np.bool_] | None = None
    if missing_data == 'flanking':
        query_valid_ranges, context_valid_mask = flanking_masks(queries, context_matrix)
    
    # Find neighbors for each query strain
    query_items = list(queries.items())
    nthreads = min(args.nthreads, len(query_items))  # don't spawn more threads than queries

    cancel_event = threading.Event()

    _process_queries_kwargs = dict(
        context_matrix=context_matrix,
        context_names=context_names,
        k=args.k,
        max_distance=args.max_distance,
        missing_data=missing_data,
        query_valid_ranges=query_valid_ranges,
        context_valid_mask=context_valid_mask,
        on_query_done=_make_progress_callback(len(query_items), not args.no_progress),
        cancel_event=cancel_event,
    )

    if nthreads <= 1:
        # Fast path: no threading overhead
        results = _process_queries(query_items, **_process_queries_kwargs)
    else:
        # Chunk queries across threads
        chunk_size = math.ceil(len(query_items) / nthreads)
        chunks = [query_items[i:i + chunk_size] for i in range(0, len(query_items), chunk_size)]
        print_err(f"Processing {len(query_items):,} queries across {len(chunks):,} threads")

        results = {}
        with ThreadPoolExecutor(max_workers=nthreads) as executor:
            futures = [executor.submit(_process_queries, chunk, **_process_queries_kwargs) for chunk in chunks]
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
    print_err(f"Found {len(all_neighbour_strain_names):,} unique neighbor strains for {len(queries)} query strains")

    # Write output in deterministic (sorted) order
    with open(args.output_strains, 'w') as fh:
        for strain_name in sorted(all_neighbour_strain_names):
            print(strain_name, file=fh)
    print_err(f"Wrote {len(all_neighbour_strain_names):,} strains to {args.output_strains}")

    # Write per-query matches TSV if requested
    if args.output_matches:
        with open(args.output_matches, 'w') as fh:
            print("query_strain\tcontext_strain\tdistance", file=fh)
            for query_strain, matches in results.items():
                for match in matches:
                    print(f"{query_strain}\t{match['strain']}\t{match['distance']}", file=fh)
        print_err(f"Wrote matches to {args.output_matches}")
