import datetime

from augur.filter import Filterer
from augur.filtering.matchers import MATCHER_CLASSES
from augur.sequence import Sequence

import argparse
import pytest


@pytest.fixture
def mock_args(mocker):
    # TODO needed?
    def _mock_args(args_dict):
        return argparse.Namespace(**args_dict)

    return _mock_args


def build_args(**args):
    return argparse.Namespace(
        **{
            "sequences": "tests/builds/tb/data/lee_2015.vcf",
            "metadata": "tests/builds/tb/data/meta.tsv",
            "min_date": "2000-05-05",
            "max_date": "2000-06-06",
            "min_length": "5",
            "non_nucleotide": True,
            "exclude": "tests/builds/tb/data/dropped_strains.txt",
            "exclude_where": "key!=val",
            "include": "tests/builds/tb/data/dropped_strains.txt",
            "include_where": "key=val",
            **args,
        }
    )


class TestFilterer:
    def test_run(self):
        # TODO end-to-end
        pass

    def test_enabled_matchers_exclude(self):
        exclude_matchers = Filterer(build_args()).enabled_matchers()["exclude"]

        assert len(exclude_matchers) == 5
        assert exclude_matchers[0].__class__ == MATCHER_CLASSES["date"]
        assert exclude_matchers[0].min_date == datetime.date(2000, 5, 5)
        assert exclude_matchers[0].max_date == datetime.date(2000, 6, 6)
        assert exclude_matchers[1].__class__ == MATCHER_CLASSES["length"]
        assert exclude_matchers[1].min_length == 5
        assert exclude_matchers[2].__class__ == MATCHER_CLASSES["non-nucleotide"]
        assert exclude_matchers[3].__class__ == MATCHER_CLASSES["name"]
        assert exclude_matchers[3].names == {"G22696"}
        assert exclude_matchers[4].__class__ == MATCHER_CLASSES["metadata"]
        assert exclude_matchers[4].conditions == [("key", "!=", "val")]

    def test_enabled_matchers_include(self):
        include_matchers = Filterer(build_args()).enabled_matchers()["include"]

        assert len(include_matchers) == 2
        assert include_matchers[0].__class__ == MATCHER_CLASSES["name"]
        assert include_matchers[0].names == {"G22696"}
        assert include_matchers[1].__class__ == MATCHER_CLASSES["metadata"]
        assert include_matchers[1].conditions == [("key", "=", "val")]

    def test_filter_sequences(self):
        pass

    def test_no_exclusion_matches(self):
        pass

    def test_any_inclusion_matches(self):
        pass

    def test_subsample(self):
        # TODO
        pass

    def test_write_output_file_vcf(self, mocker, sequence_factory):
        samples = sequence_factory.build_batch(3)

        filterer = Filterer(build_args(output="output.vcf"))
        filterer.resulting_sequences = samples
        filterer.sequence_file = mocker.Mock()
        vcf_selector_class = mocker.patch("augur.filter.VCFSelector")

        filterer.write_output_file()

        vcf_selector_class.assert_called_once_with(
            selected_samples=samples,
            samplefile=filterer.sequence_file,
            output_fname="output.vcf",
        )

    def test_write_output_file_fasta(self, mocker, sequence_factory):
        samples = sequence_factory.build_batch(3)
        filterer = Filterer(build_args(sequences="tests/builds/tb/data/ref.fasta", output="output.fasta"))
        filterer.resulting_sequences = samples
        bio_seqio = mocker.patch("augur.filter.SeqIO")

        filterer.write_output_file()

        bio_seqio.write.assert_called_once_with(
            [sample.bio_seq for sample in samples],
            "output.fasta",
            "fasta",
        )

    def test_print_summary(self):
        pass

    def check_vcftools(self):
        pass
