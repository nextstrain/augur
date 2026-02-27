"""Tests for _process_queries"""
import numpy as np

from augur.proximity import (
    _process_queries,
    get_valid_range,
    to_numpy_array,
)


def _build(seqs: dict[str, str]):
    """Helper: convert string sequences to context_matrix + context_names."""
    names = list(seqs.keys())
    matrix = np.stack([to_numpy_array(s) for s in seqs.values()])
    return matrix, names


def _get_neighbours(query_seq, context_matrix, context_names, k, max_distance,
                    missing_data, query_valid_range=None, context_valid_mask=None):
    """Thin wrapper: run _process_queries for a single query named '_q'."""
    query_valid_ranges = {"_q": query_valid_range} if query_valid_range else None
    results = _process_queries(
        [("_q", query_seq)],
        context_matrix, context_names, k, max_distance, missing_data,
        query_valid_ranges=query_valid_ranges,
        context_valid_mask=context_valid_mask,
    )
    return results["_q"]


class TestMissingDataNone:
    """missing_data='none' — N is treated as a regular base."""

    def test_basic_distances(self):
        query = to_numpy_array("ATCGATCG")
        context_matrix, context_names = _build({
            "exact":  "ATCGATCG",
            "one":    "TTCGATCG",  # pos 0 differs
            "two":    "TTCGTTCG",  # pos 0,4 differ
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=10, max_distance=10, missing_data="none",
        )
        assert result == [{"strain": "exact", "distance": 0}, {"strain": "one", "distance": 1}, {"strain": "two", "distance": 2}]

    def test_n_counted_as_mismatch(self):
        query = to_numpy_array("ATCG")
        context_matrix, context_names = _build({
            "has_n": "ANCG",  # N at pos 1 counts as mismatch
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=10, max_distance=10, missing_data="none",
        )
        assert result == [{"strain": "has_n", "distance": 1}]

    def test_k_limits_results(self):
        query = to_numpy_array("ATCGATCG")
        context_matrix, context_names = _build({
            "d0": "ATCGATCG",
            "d1": "TTCGATCG",
            "d2": "TTCGTTCG",
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=2, max_distance=10, missing_data="none",
        )
        assert result == [{"strain": "d0", "distance": 0}, {"strain": "d1", "distance": 1}]

    def test_max_distance_filters(self):
        query = to_numpy_array("ATCGATCG")
        context_matrix, context_names = _build({
            "d1": "TTCGATCG",
            "d3": "TTCGTTCG",  # 2 mismatches, outside threshold
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=10, max_distance=1, missing_data="none",
        )
        assert result == [{"strain": "d1", "distance": 1}]

    def test_no_matches(self):
        query = to_numpy_array("ATCG")
        context_matrix, context_names = _build({
            "far": "TAGC",  # 4 mismatches
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=5, max_distance=1, missing_data="none",
        )
        assert result == []

    def test_tiebreak_fewest_ns(self):
        """When distances are equal, prefer sequences with fewer Ns."""
        query = to_numpy_array("ATCGATCG")
        context_matrix, context_names = _build({
            "many_n":  "TTCGANNN",  # distance 4 (pos 0 + 3 Ns), 3 Ns
            "few_n":   "TTCGATNG",  # distance 2 (pos 0 + 1 N), 1 N
            "no_n":    "TTCGATCG",  # distance 1 (pos 0), 0 Ns
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=3, max_distance=10, missing_data="none",
        )
        # sorted by (distance, n_count, name): no_n(1,0), few_n(2,1), many_n(4,3)
        assert result == [{"strain": "no_n", "distance": 1}, {"strain": "few_n", "distance": 2}, {"strain": "many_n", "distance": 4}]

    def test_tiebreak_fewest_ns_same_distance(self):
        """When distances are equal, fewer Ns wins over alphabetical order."""
        query = to_numpy_array("ATCGATCG")
        context_matrix, context_names = _build({
            # All have distance 2 from query (pos 0 mismatch + N mismatches)
            "alpha":   "TTCNANCG",  # distance 3, 2 Ns
            "bravo":   "TTCGATNG",  # distance 2, 1 N
            "charlie": "TTCGATCG",  # distance 1, 0 Ns
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=3, max_distance=10, missing_data="none",
        )
        # charlie(1, 0 Ns) < bravo(2, 1 N) < alpha(3, 2 Ns)
        assert result == [{"strain": "charlie", "distance": 1}, {"strain": "bravo", "distance": 2}, {"strain": "alpha", "distance": 3}]

    def test_tiebreak_fewest_ns_at_k_boundary(self):
        """At the k boundary with equal distances, fewer Ns wins over alphabetical."""
        query = to_numpy_array("ATCGATCG")
        #                       01234567
        context_matrix, context_names = _build({
            "best":    "TTCGATCG",  # distance 1, 0 Ns
            "clean":   "TTCGTTCG",  # distance 2, 0 Ns
            "has_n":   "TNCGATCG",  # distance 2, 1 N  (pos 0: A!=T, pos 1: T!=N)
            "more_n":  "NNCGATCG",  # distance 2, 2 Ns (pos 0: A!=N, pos 1: T!=N)
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=3, max_distance=10, missing_data="none",
        )
        # best(dist=1) is below cutoff. Cutoff distance is 2 with 3 candidates, need 2.
        # At cutoff, sorted by (n_count, name): clean(0), has_n(1), more_n(2) -> take clean, has_n
        assert result == [{"strain": "best", "distance": 1}, {"strain": "clean", "distance": 2}, {"strain": "has_n", "distance": 2}]

    def test_tiebreak_alphabetical_fallback(self):
        """When distances and N counts are equal, fall back to alphabetical."""
        query = to_numpy_array("ATCG")
        context_matrix, context_names = _build({
            "zulu":  "TTCG",  # distance 1, 0 Ns
            "alpha": "TTCG",  # distance 1, 0 Ns
            "mike":  "TTCG",  # distance 1, 0 Ns
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=2, max_distance=10, missing_data="none",
        )
        assert result == [{"strain": "alpha", "distance": 1}, {"strain": "mike", "distance": 1}]


class TestMissingDataAll:
    """missing_data='all' — positions where either sequence has N are ignored."""

    def test_query_n_ignored(self):
        query = to_numpy_array("ANCG")  # N at pos 1
        context_matrix, context_names = _build({
            "seq": "ATCG",
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=10, max_distance=10, missing_data="all",
        )
        assert result == [{"strain": "seq", "distance": 0}]

    def test_context_n_ignored(self):
        query = to_numpy_array("ATCG")
        context_matrix, context_names = _build({
            "seq": "ANCG",  # N at pos 1
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=10, max_distance=10, missing_data="all",
        )
        assert result == [{"strain": "seq", "distance": 0}]

    def test_both_n_ignored(self):
        query = to_numpy_array("ANNG")
        context_matrix, context_names = _build({
            "seq": "ANCG",  # both have N at pos 1; context N at pos 2
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=10, max_distance=10, missing_data="all",
        )
        # pos 1: both N -> ignored; pos 2: query N -> ignored
        assert result == [{"strain": "seq", "distance": 0}]

    def test_mismatch_on_valid_positions(self):
        query = to_numpy_array("ATCG")
        context_matrix, context_names = _build({
            "one_real": "ANCG",  # pos 1 is N (ignored), pos 0 matches -> distance 0
            "one_mut":  "TTCG",  # pos 0 differs -> distance 1
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=10, max_distance=10, missing_data="all",
        )
        assert result == [{"strain": "one_real", "distance": 0}, {"strain": "one_mut", "distance": 1}]

    def test_all_n_gives_zero_distance(self):
        query = to_numpy_array("NNNN")
        context_matrix, context_names = _build({
            "seq": "ATCG",
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=10, max_distance=10, missing_data="all",
        )
        assert result == [{"strain": "seq", "distance": 0}]


class TestMissingDataFlanking:
    """missing_data='flanking' — runs of N at start/end of each sequence are ignored."""

    @staticmethod
    def _build_flanking(context_seqs: dict[str, str]):
        """Build context_matrix and context_valid_mask for flanking mode."""
        names = list(context_seqs.keys())
        arrays = [to_numpy_array(s) for s in context_seqs.values()]
        matrix = np.stack(arrays)
        mask = np.ones_like(matrix, dtype=bool)
        for i, arr in enumerate(arrays):
            start, end = get_valid_range(arr)
            mask[i, :start] = False
            mask[i, end:] = False
        return matrix, names, mask

    def test_flanking_ns_ignored(self):
        query = to_numpy_array("ATCGATCG")
        q_range = get_valid_range(query)
        context_matrix, context_names, mask = self._build_flanking({
            "seq": "NNATCGNN",  # flanking Ns at pos 0-1 and 6-7
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=10, max_distance=10,
            missing_data="flanking", query_valid_range=q_range, context_valid_mask=mask,
        )
        # Only compare pos 2-5 (CGAT vs CGAT in context after lowering) -> distance 0
        # Wait, let me think: query="ATCGATCG", context="NNATCGNN"
        # context valid range is [2, 6), query valid range is [0, 8)
        # intersection is [2, 6), comparing query[2:6]="CGAT" vs context[2:6]="ATCG"
        # C!=A, G!=T, A!=C, T!=G -> 4 mismatches
        assert result == [{"strain": "seq", "distance": 4}]

    def test_flanking_ns_both_sides(self):
        query = to_numpy_array("NNATCGNN")
        q_range = get_valid_range(query)
        context_matrix, context_names, mask = self._build_flanking({
            "seq": "NNATCGNN",
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=10, max_distance=10,
            missing_data="flanking", query_valid_range=q_range, context_valid_mask=mask,
        )
        # Both have valid range [2, 6), same content -> distance 0
        assert result == [{"strain": "seq", "distance": 0}]

    def test_overlapping_valid_ranges(self):
        query = to_numpy_array("NATCGATN")  # valid [1, 7)
        q_range = get_valid_range(query)
        context_matrix, context_names, mask = self._build_flanking({
            "seq": "NNNTCGNN",  # valid [3, 6)
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=10, max_distance=10,
            missing_data="flanking", query_valid_range=q_range, context_valid_mask=mask,
        )
        # intersection is [3, 6), query[3:6]="GAT", context[3:6]="TCG"
        # G!=T, A!=C, T!=G -> 3 mismatches
        assert result == [{"strain": "seq", "distance": 3}]

    def test_no_overlap_excluded(self):
        query = to_numpy_array("ATCNNNNN")  # valid [0, 3)
        q_range = get_valid_range(query)
        context_matrix, context_names, mask = self._build_flanking({
            "seq": "NNNNATCG",  # valid [4, 8)
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=10, max_distance=10,
            missing_data="flanking", query_valid_range=q_range, context_valid_mask=mask,
        )
        # No overlap -> 0 mismatches on 0 positions -> distance 0
        assert result == [{"strain": "seq", "distance": 0}]

    def test_interior_ns_still_counted(self):
        """Flanking mode only ignores leading/trailing Ns, not interior ones."""
        query = to_numpy_array("ATCGATCG")
        q_range = get_valid_range(query)
        context_matrix, context_names, mask = self._build_flanking({
            "seq": "ATNGNTCG",  # interior Ns at pos 2, 4
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=10, max_distance=10,
            missing_data="flanking", query_valid_range=q_range, context_valid_mask=mask,
        )
        # context valid range [0, 8), query valid range [0, 8)
        # pos 2: C!=N, pos 4: A!=N -> 2 mismatches
        assert result == [{"strain": "seq", "distance": 2}]

    def test_query_flanking_ns(self):
        """Query's flanking Ns should also narrow the comparison range."""
        query = to_numpy_array("NNATCGNN")  # valid [2, 6)
        q_range = get_valid_range(query)
        context_matrix, context_names, mask = self._build_flanking({
            "exact": "TTATCGCC",  # valid [0, 8), but only [2,6) compared
            "diff":  "TTTTCGCC",  # valid [0, 8), pos 2: A!=T -> 1 mismatch in overlap
        })
        result = _get_neighbours(
            query, context_matrix, context_names, k=10, max_distance=10,
            missing_data="flanking", query_valid_range=q_range, context_valid_mask=mask,
        )
        assert result == [{"strain": "exact", "distance": 0}, {"strain": "diff", "distance": 1}]
