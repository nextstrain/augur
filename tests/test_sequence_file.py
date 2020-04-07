from augur.sequence_file import SequenceFile

import pytest


@pytest.fixture
def mock_metadata(mocker):
    def _mock_metadata(names):
        mocker.patch(
            "augur.sequence_file.read_metadata",
            lambda x: ({name: {"name": name, "date": "5"} for name in names}, []),
        )

    return _mock_metadata


class TestSequenceFile:
    def test_sequences_vcf(self, mock_metadata):
        mock_metadata(["G22565", "G22566", "G22567", "G22568", "G22569", "G22571"])

        sequences = SequenceFile("tests/builds/tb/data/lee_2015.vcf", "").sequences

        assert len(sequences) == 6
        assert sequences[5].name == "G22571"
        assert sequences[5].metadata["date"] == "5"

    def test_sequences_vcf_gz(self, mock_metadata):
        mock_metadata(["G22565", "G22566", "G22567", "G22568", "G22569", "G22571"])

        sequences = SequenceFile("tests/builds/tb/data/lee_2015.vcf.gz", "").sequences

        assert len(sequences) == 6
        assert sequences[5].name == "G22571"
        assert sequences[5].metadata["date"] == "5"

    def test_sequences_fasta(self, mock_metadata):
        mock_metadata(["AAK51718.1"])

        sequences = SequenceFile(
            "tests/data/fitness_model/AAK51718.fasta", ""
        ).sequences

        assert len(sequences) == 1
        assert str(sequences[0].sequence)[:5] == "MKTII"
        assert sequences[0].metadata["date"] == "5"

    def test_sequences_unsupported_format(self, mocker):
        # TODO make classes for all of these exceptions
        with pytest.raises(Exception):
            SequenceFile("path/to/file.json", "").sequences
