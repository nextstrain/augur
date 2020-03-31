import augur.filter
import pytest


@pytest.fixture
def mock_priorities_file_valid(mocker):
    mocker.patch(
        "builtins.open", mocker.mock_open(read_data="strain1 5\nstrain2 6\nstrain3 8\n")
    )


@pytest.fixture
def mock_priorities_file_malformed(mocker):
    mocker.patch("builtins.open", mocker.mock_open(read_data="strain1 X\n"))


class TestFilter:
    def test_read_vcf_compressed(self):
        seq_keep, all_seq = augur.filter.read_vcf(
            "tests/builds/tb/data/lee_2015.vcf.gz"
        )

        assert len(seq_keep) == 150
        assert seq_keep[149] == "G22733"
        assert seq_keep == all_seq

    def test_read_vcf_uncompressed(self):
        seq_keep, all_seq = augur.filter.read_vcf("tests/builds/tb/data/lee_2015.vcf")

        assert len(seq_keep) == 150
        assert seq_keep[149] == "G22733"
        assert seq_keep == all_seq

    def test_read_priority_scores_valid(self, mock_priorities_file_valid):
        # builtins.open is stubbed, but we need a valid file to satisfy the existence check
        priorities = augur.filter.read_priority_scores(
            "tests/builds/tb/data/lee_2015.vcf"
        )

        assert priorities == {"strain1": 5, "strain2": 6, "strain3": 8}

    def test_read_priority_scores_malformed(self, mock_priorities_file_malformed):
        with pytest.raises(ValueError):
            # builtins.open is stubbed, but we need a valid file to satisfy the existence check
            augur.filter.read_priority_scores("tests/builds/tb/data/lee_2015.vcf")

    def test_read_priority_scores_does_not_exist(self):
        with pytest.raises(FileNotFoundError):
            augur.filter.read_priority_scores("/does/not/exist.txt")
