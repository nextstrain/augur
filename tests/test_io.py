#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import bz2
import gzip
import lzma
from pathlib import Path
import pytest
import random
import sys

import augur.io as io
from augur.io import open_file, read_sequences, write_sequences


def random_seq(k):
    """Generate a single random sequence of nucleotides of length k.
    """
    return "".join(random.choices(("A","T","G","C"), k=k))

def generate_sequences(n, k=10):
    """Generate n random sequences of length k.
    """
    return (
        SeqRecord(Seq(random_seq(k)), id=f"SEQ_{i}")
        for i in range(1, n + 1)
    )

@pytest.fixture
def sequences():
    return list(generate_sequences(3))

@pytest.fixture
def sequences_generator():
    return generate_sequences(3)

@pytest.fixture
def fasta_filename(tmpdir, sequences):
    filename = str(tmpdir / "sequences.fasta")
    SeqIO.write(sequences, filename, "fasta")
    return filename

@pytest.fixture
def additional_fasta_filename(tmpdir, sequences):
    filename = str(tmpdir / "additional_sequences.fasta")
    SeqIO.write(sequences, filename, "fasta")
    return filename

@pytest.fixture
def gzip_fasta_filename(tmpdir, sequences):
    filename = str(tmpdir / "sequences.fasta.gz")

    with gzip.open(filename, "wt") as oh:
        SeqIO.write(sequences, oh, "fasta")

    return filename

@pytest.fixture
def bzip2_fasta_filename(tmpdir, sequences):
    filename = str(tmpdir / "sequences.fasta.bz2")

    with bz2.open(filename, "wt") as oh:
        SeqIO.write(sequences, oh, "fasta")

    return filename

@pytest.fixture
def lzma_fasta_filename(tmpdir, sequences):
    filename = str(tmpdir / "sequences.fasta.xz")

    with lzma.open(filename, "wt") as oh:
        SeqIO.write(sequences, oh, "fasta")

    return filename

@pytest.fixture
def genbank_reference():
    return "tests/builds/zika/config/zika_outgroup.gb"

@pytest.fixture
def mock_run_shell_command(mocker):
    mocker.patch("augur.io.run_shell_command")


class TestReadSequences:
    def test_read_sequences_from_single_file(self, fasta_filename):
        sequences = read_sequences(fasta_filename, format="fasta")
        assert len(list(sequences)) == 3

    def test_read_sequences_from_multiple_files(self, fasta_filename, additional_fasta_filename):
        sequences = read_sequences(fasta_filename, additional_fasta_filename, format="fasta")
        assert len(list(sequences)) == 6

    def test_read_sequences_from_multiple_files_or_buffers(self, fasta_filename, additional_fasta_filename):
        with open(fasta_filename) as fasta_handle:
            sequences = read_sequences(fasta_handle, additional_fasta_filename, format="fasta")
            assert len(list(sequences)) == 6

    def test_read_single_fasta_record(self, fasta_filename):
        record = next(read_sequences(fasta_filename, format="fasta"))
        assert record.id == "SEQ_1"

    def test_read_single_genbank_record(self, genbank_reference):
        reference = next(read_sequences(genbank_reference, format="genbank"))
        assert reference.id == "KX369547.1"

    def test_read_single_genbank_record_from_a_path(self, genbank_reference):
        reference = next(read_sequences(Path(genbank_reference), format="genbank"))
        assert reference.id == "KX369547.1"

    def test_read_sequences_from_single_gzip_file(self, gzip_fasta_filename):
        sequences = read_sequences(gzip_fasta_filename, format="fasta")
        assert len(list(sequences)) == 3

    def test_read_sequences_from_single_lzma_file(self, lzma_fasta_filename):
        sequences = read_sequences(lzma_fasta_filename, format="fasta")
        assert len(list(sequences)) == 3

    def test_read_sequences_from_single_bzip2_file(self, bzip2_fasta_filename):
        sequences = read_sequences(bzip2_fasta_filename, format="fasta")
        assert len(list(sequences)) == 3

    def test_read_sequences_from_multiple_files_with_different_compression(self, fasta_filename, gzip_fasta_filename, lzma_fasta_filename):
        sequences = read_sequences(fasta_filename, gzip_fasta_filename, lzma_fasta_filename, format="fasta")
        assert len(list(sequences)) == 9


class TestWriteSequences:
    def test_write_sequences(self, tmpdir, sequences):
        output_filename = Path(tmpdir) / Path("new_sequences.fasta")
        sequences_written = write_sequences(sequences, output_filename, "fasta")
        assert sequences_written == len(sequences)

    def test_write_genbank_sequence(self, tmpdir, genbank_reference):
        output_filename = Path(tmpdir) / Path("new_sequences.fasta")

        reference = SeqIO.read(genbank_reference, "genbank")
        sequences_written = write_sequences([reference], output_filename, "genbank")
        assert sequences_written == 1

    def test_write_sequences_from_generator(self, tmpdir, sequences_generator):
        output_filename = Path(tmpdir) / Path("new_sequences.fasta")
        sequences_written = write_sequences(sequences_generator, output_filename, "fasta")
        assert sequences_written == 3

    def test_write_single_set_of_sequences_to_gzip_file(self, tmpdir, sequences):
        output_filename = Path(tmpdir) / Path("new_sequences.fasta.gz")
        sequences_written = write_sequences(sequences, output_filename, "fasta")
        assert sequences_written == len(sequences)

        with gzip.open(output_filename, "rt") as handle:
            assert sequences_written == len([line for line in handle if line.startswith(">")])

    def test_write_single_set_of_sequences_to_bzip2_file(self, tmpdir, sequences):
        output_filename = Path(tmpdir) / Path("new_sequences.fasta.bz2")
        sequences_written = write_sequences(sequences, output_filename, "fasta")
        assert sequences_written == len(sequences)

        with bz2.open(output_filename, "rt") as handle:
            assert sequences_written == len([line for line in handle if line.startswith(">")])

    def test_write_single_set_of_sequences_to_lzma_file(self, tmpdir, sequences):
        output_filename = Path(tmpdir) / Path("new_sequences.fasta.xz")
        sequences_written = write_sequences(sequences, output_filename, "fasta")
        assert sequences_written == len(sequences)

        with lzma.open(output_filename, "rt") as handle:
            assert sequences_written == len([line for line in handle if line.startswith(">")])

    def test_write_sequences_by_external_handle(self, tmpdir, sequences):
        output_filename = Path(tmpdir) / Path("new_sequences.fasta")

        with open_file(output_filename, "w") as handle:
            total_sequences_written = 0
            for sequence in sequences:
                sequences_written = write_sequences(
                    sequence,
                    handle
                )
                total_sequences_written += sequences_written

        with open(output_filename, "r") as handle:
            assert total_sequences_written == len([line for line in handle if line.startswith(">")])


class TestOpenFile:
    def test_open_file_read_text(self, tmpdir):
        """Read a text file."""
        path = str(tmpdir / 'test.txt')
        with open(path, 'w') as f_write:
            f_write.write('foo\nbar\n')
        with open_file(path) as f_read:
            assert f_read.read() == 'foo\nbar\n'

    def test_open_file_write_text(self, tmpdir):
        """Write a text file."""
        path = str(tmpdir / 'test.txt')
        with open_file(path, 'w') as f_write:
            f_write.write('foo\nbar\n')
        with open(path) as f_read:
            assert f_read.read() == 'foo\nbar\n'

    def test_open_file_read_gzip(self, tmpdir):
        """Read a text file compressed with gzip."""
        import gzip
        path = str(tmpdir / 'test.txt.gz')
        with gzip.open(path, 'wt') as f_write:
            f_write.write('foo\nbar\n')
        with open_file(path) as f_read:
            assert f_read.read() == 'foo\nbar\n'

    def test_open_file_write_gzip(self, tmpdir):
        """Write a text file compressed with gzip."""
        import gzip
        path = str(tmpdir / 'test.txt.gz')
        with open_file(path, 'w') as f_write:
            f_write.write('foo\nbar\n')
        with gzip.open(path, 'rt') as f_read:
            assert f_read.read() == 'foo\nbar\n'

    def test_open_file_read_lzma(self, tmpdir):
        """Read a text file compressed with LZMA."""
        import lzma
        path = str(tmpdir / 'test.txt.xz')
        with lzma.open(path, 'wt') as f_write:
            f_write.write('foo\nbar\n')
        with open_file(path) as f_read:
            assert f_read.read() == 'foo\nbar\n'

    def test_open_file_read_lzma(self, tmpdir):
        """Write a text file compressed with LZMA."""
        import lzma
        path = str(tmpdir / 'test.txt.xz')
        with open_file(path, 'w') as f_write:
            f_write.write('foo\nbar\n')
        with lzma.open(path, 'rt') as f_read:
            assert f_read.read() == 'foo\nbar\n'


class TestVCF:
    def test_read_vcf_compressed(self):
        seq_keep, all_seq = io.read_vcf(
            "tests/builds/tb/data/lee_2015.vcf.gz"
        )

        assert len(seq_keep) == 150
        assert seq_keep[149] == "G22733"
        assert seq_keep == all_seq

    def test_read_vcf_uncompressed(self):
        seq_keep, all_seq = io.read_vcf("tests/builds/tb/data/lee_2015.vcf")

        assert len(seq_keep) == 150
        assert seq_keep[149] == "G22733"
        assert seq_keep == all_seq

    def test_write_vcf_compressed_input(self, mock_run_shell_command):
        io.write_vcf(
            "tests/builds/tb/data/lee_2015.vcf.gz", "output_file.vcf.gz", []
        )

        io.run_shell_command.assert_called_once_with(
            "vcftools --gzvcf tests/builds/tb/data/lee_2015.vcf.gz --recode --stdout | gzip -c > output_file.vcf.gz",
            raise_errors=True,
        )

    def test_write_vcf_uncompressed_input(self, mock_run_shell_command):
        io.write_vcf(
            "tests/builds/tb/data/lee_2015.vcf", "output_file.vcf.gz", []
        )

        io.run_shell_command.assert_called_once_with(
            "vcftools --vcf tests/builds/tb/data/lee_2015.vcf --recode --stdout | gzip -c > output_file.vcf.gz",
            raise_errors=True,
        )

    def test_write_vcf_compressed_output(self, mock_run_shell_command):
        io.write_vcf(
            "tests/builds/tb/data/lee_2015.vcf", "output_file.vcf.gz", []
        )

        io.run_shell_command.assert_called_once_with(
            "vcftools --vcf tests/builds/tb/data/lee_2015.vcf --recode --stdout | gzip -c > output_file.vcf.gz",
            raise_errors=True,
        )

    def test_write_vcf_uncompressed_output(self, mock_run_shell_command):
        io.write_vcf(
            "tests/builds/tb/data/lee_2015.vcf", "output_file.vcf", []
        )

        io.run_shell_command.assert_called_once_with(
            "vcftools --vcf tests/builds/tb/data/lee_2015.vcf --recode --stdout  > output_file.vcf",
            raise_errors=True,
        )

    def test_write_vcf_dropped_samples(self, mock_run_shell_command):
        io.write_vcf(
            "tests/builds/tb/data/lee_2015.vcf", "output_file.vcf", ["x", "y", "z"]
        )

        io.run_shell_command.assert_called_once_with(
            "vcftools --remove-indv x --remove-indv y --remove-indv z --vcf tests/builds/tb/data/lee_2015.vcf --recode --stdout  > output_file.vcf",
            raise_errors=True,
        )
