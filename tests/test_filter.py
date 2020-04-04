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
            "min_date": None,
            "max_date": None,
            "min_length": None,
            "non_nucleotide": None,
            "exclude": None,
            "exclude_where": None,
            "include": None,
            "include_where": None,
            **args
        }
    )


@pytest.fixture
def mock_translate_args(mocker):
    mocker.patch("augur.filter.Filterer.translate_args", lambda _: {
        "exclude": ["metadata:key!=val", "name:file=tests/builds/tb/data/dropped_strains.txt"],
        "include": ["length:min=5"],
    })


class TestFilterer:
    def test_run(self):
        # TODO end-to-end
        pass

    @pytest.mark.parametrize(
        "args, expected_exclude, expected_include",
        [
            (
                build_args(min_length=5, min_date=2000),
                ["date:min=2000", "length:min=5"],
                [],
            ),
            (
                build_args(non_nucleotide=True, include="include_file.txt"),
                ["non-nucleotide"],
                ["name:file=include_file.txt"],
            ),
            (
                build_args(exclude="exclude_file.txt", include_where="key=val"),
                ["name:file=exclude_file.txt"],
                ["metadata:key=val"],
            ),
            (
                build_args(exclude_where="key!=val"),
                ["metadata:key!=val"],
                [],
            ),
        ]
    )
    def test_translate_args(self, args, expected_exclude, expected_include):
        assert Filterer(args).translate_args() == {
            "exclude": expected_exclude,
            "include": expected_include,
        }

    def test_enabled_matchers_exclude(self, mock_translate_args):
        enabled_matchers = Filterer(build_args()).enabled_matchers("exclude")

        assert len(enabled_matchers) == 2
        assert enabled_matchers[0].__class__ == MATCHER_CLASSES["metadata"]
        assert enabled_matchers[0].conditions == [("key", "!=", "val")]
        assert enabled_matchers[1].__class__ == MATCHER_CLASSES["name"]
        assert enabled_matchers[1].names == {"G22696"}

    def test_enabled_matchers_include(self, mock_translate_args):
        enabled_matchers = Filterer(build_args()).enabled_matchers("include")

        assert len(enabled_matchers) == 1
        assert enabled_matchers[0].__class__ == MATCHER_CLASSES["length"]
        assert enabled_matchers[0].min_length == 5

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
