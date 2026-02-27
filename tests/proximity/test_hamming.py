"""Tests for augur proximity (hamming distance method)"""
from __future__ import annotations
import argparse
import lzma
import pytest

from augur import proximity

class Genome:
    def __init__(self, starting_genome: str) -> None:
        self.seq = list(starting_genome)
        self.mutated_positions: set[int] = set()
        self.num_ns: int = 0
    
    def __str__(self) -> str:
        return "".join(self.seq)

    def mutate(self, n_mutations: int) -> Genome:
        """Introduce *n_mutations* deterministic mutations into *seq*.
        Mutations are placed at positions 0, 4, 8, ... (every 4th base)
        Each base is swapped: A<->T, C<->G.
        """
        swap = {"A": "T", "T": "A", "C": "G", "G": "C"}
        for i in range(n_mutations):
            pos = i * 4
            self.seq[pos] = swap[self.seq[pos]] # raises IndexError if we mess up
            self.mutated_positions.add(pos)
        return self

    def missing(self, /, starting:int=0, trailing:int=0, num:int=0) -> Genome:
        for i in range(0, starting):
            self.seq[i] = 'N'
            self.mutated_positions.discard(i)
            self.num_ns += 1
        for i in range(len(self.seq)-trailing, len(self.seq)):
            self.seq[i] = 'N'
            self.mutated_positions.discard(i)
            self.num_ns += 1
        for i in range(num):
            pos = i*7
            if self.seq[pos]!='N':
                self.seq[pos] = 'N'
                self.mutated_positions.discard(pos)
                self.num_ns += 1
        return self

def write_fasta(path: str, records: dict[str, str]):
    with lzma.open(path, "wt") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n{seq}\n")


@pytest.fixture
def query_fasta(tmp_path):
    """Factory fixture: call with optional records dict, defaults to a single
    query sequence using BASE_GENOME."""
    def _make(records: dict[str, str]):
        path = str(tmp_path / "query.fasta.xz")
        write_fasta(path, records)
        return path
    return _make


@pytest.fixture
def sequences_fasta(tmp_path):
    """Factory fixture: call with optional records dict, defaults to 10 context
    sequences with 1-10 mutations each."""
    def _make(records: dict[str, str]):            
        path = str(tmp_path / "sequences.fasta.xz")
        write_fasta(path, records)
        return path
    return _make


@pytest.fixture
def output_file(tmp_path):
    return str(tmp_path / "output_strains.txt")


@pytest.fixture
def argparser():
    """Provide an easy way to test command line arguments."""
    parser = argparse.ArgumentParser()
    proximity.register_arguments(parser)
    def parse(args):
        return parser.parse_args(args.split())
    return parse


def read_output(path) -> list[str]:
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


class TestHammingSimple:
    BASE_GENOME = "ATCG" * 10  # 40nt
    QUERY = {"query": BASE_GENOME}
    SEQS = {}
    for i in range(1,11): # no dict comprehension - BASE_GENOME out of scope
        SEQS[f"seq_{i}"] = str(Genome(BASE_GENOME).mutate(i))

    def test_run_finds_correct_neighbours(self, query_fasta, sequences_fasta, output_file, argparser):
        """With k=5 and a high max-distance, the 5 closest sequences should be returned."""
        args = argparser(f"--method hamming --query {query_fasta(self.QUERY)} --sequences {sequences_fasta(self.SEQS)} --output-strains {output_file} --k 5 --max-distance 10")
        proximity.run(args)
        result = read_output(output_file)
        assert sorted(result) == ["seq_1", "seq_2", "seq_3", "seq_4", "seq_5"]

    def test_run_k_limits_output(self, query_fasta, sequences_fasta, output_file, argparser):
        """k=3 should return only the 3 closest sequences."""
        args = argparser(f"--method hamming --query {query_fasta(self.QUERY)} --sequences {sequences_fasta(self.SEQS)} --output-strains {output_file} --k 3 --max-distance 10")
        proximity.run(args)
        result = read_output(output_file)
        assert sorted(result) == ["seq_1", "seq_2", "seq_3"]

    def test_run_max_distance_arg(self, query_fasta, sequences_fasta, output_file, argparser):
        """With max-distance of 3 we should only get the first three sequences"""
        args = argparser(f"--method hamming --query {query_fasta(self.QUERY)} --sequences {sequences_fasta(self.SEQS)} --output-strains {output_file} --k 10 --max-distance 3")
        proximity.run(args)
        result = read_output(output_file)
        assert sorted(result) == ["seq_1", "seq_2", "seq_3"]

    def test_run_excludes_query_from_context(self, query_fasta, sequences_fasta, output_file, argparser):
        """If the query strain also appears in the context FASTA, it should be excluded."""
        # Create a context FASTA that includes the query sequence
        records = {**self.QUERY}
        records.update({**self.SEQS})
        args = argparser(f"--method hamming --query {query_fasta(self.QUERY)} --sequences {sequences_fasta(records)} --output-strains {output_file} --k 5 --max-distance 10")
        proximity.run(args)
        result = read_output(output_file)
        assert "query" not in result


class TestHammingWithMisingData:
    BASE_GENOME = "ATCG" * 10  # 40nt
    QUERY = {"query": BASE_GENOME}
    SEQS = {}
    SEQS_INFO = {}
    for i in range(1,6):
        SEQS_INFO[f"seq_{i}"] = Genome(BASE_GENOME).mutate(i).missing(starting=3, trailing=3, num=2)
        SEQS[f"seq_{i}"] = str(SEQS_INFO[f"seq_{i}"])
        
    def test_count_all_ns(self, query_fasta, sequences_fasta, output_file, argparser):
        max_dist = 10
        args = argparser(f"--method hamming --missing-data none --query {query_fasta(self.QUERY)} --sequences {sequences_fasta(self.SEQS)} --output-strains {output_file} --k 100 --max-distance {max_dist}")
        expected = [k for k,v in self.SEQS_INFO.items() if len(v.mutated_positions)+v.num_ns <= max_dist]
        proximity.run(args)
        result = read_output(output_file)
        assert sorted(result) == sorted(expected)

    def test_ignore_all_ns(self, query_fasta, sequences_fasta, output_file, argparser):
        max_dist = 3
        args = argparser(f"--method hamming --missing-data all --query {query_fasta(self.QUERY)} --sequences {sequences_fasta(self.SEQS)} --output-strains {output_file} --k 100 --max-distance {max_dist}")
        expected = [k for k,v in self.SEQS_INFO.items() if len(v.mutated_positions) <= max_dist]
        proximity.run(args)
        result = read_output(output_file)
        assert sorted(result) == sorted(expected)

    def test_ignore_flanking_ns(self, query_fasta, sequences_fasta, output_file, argparser):
        max_dist = 3
        args = argparser(f"--method hamming --missing-data flanking --query {query_fasta(self.QUERY)} --sequences {sequences_fasta(self.SEQS)} --output-strains {output_file} --k 100 --max-distance {max_dist}")
        # 6 flanking Ns which we exclude here
        expected = [k for k,v in self.SEQS_INFO.items() if len(v.mutated_positions)+v.num_ns-6 <= max_dist]
        proximity.run(args)
        result = read_output(output_file)
        assert sorted(result) == sorted(expected)
