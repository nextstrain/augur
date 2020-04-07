"""
Filter and subsample a sequence set.
"""
from Bio import SeqIO
import functools

from augur.exceptions import MissingDependencyException
from augur.filtering.matchers import MATCHER_CLASSES
from augur.filtering.vcf_selector import VCFSelector
from augur.sequence_file import SequenceFile


def register_arguments(parser):
    parser.add_argument(
        "--sequences",
        "-s",
        required=True,
        help="sequences in FASTA, VCF, or gzipped VCF format (detected by file extension)",
    )
    parser.add_argument(
        "--metadata", required=True, help="metadata associated with sequences"
    )
    parser.add_argument("--output", "-o", required=True, help="output file")

    parser.add_argument(
        "--sequences-per-group",
        type=int,
        help="subsample to no more than this number of sequences per category",
    )
    parser.add_argument(
        "--group-by",
        nargs="+",
        help="categories with respect to subsample; two virtual fields, `month` and `year`, are supported if they don't already exist as real fields but a `date` field does exist",
    )
    parser.add_argument(
        "--priority",
        type=str,
        help="file with list priority scores for sequences (strain\tpriority)",
    )
    parser.add_argument(
        "--subsample-seed",
        help="random number generator seed to allow reproducible sub-sampling (with same input data). Can be number or string.",
    )

    for cls in MATCHER_CLASSES.values():
        cls.add_arguments(parser)


def run(args):
    """
    Apply exclusion, inclusion, and subsampling rules to produce a subset of input data.

    This command supports exclusion or inclusion rules TODO TODO TODO

    Also subsampling TODO
    """
    Filterer(args).run()


class Filterer:
    def __init__(self, args):
        self.args = args
        self.sequence_file = SequenceFile(args.sequences, args.metadata)

        self.check_vcftools()

        self.sequences = self.sequence_file.sequences

    def run(self):
        self.filter_sequences()
        self.print_summary()
        self.write_output_file()

    def translate_args(self):
        """
        Translate legacy CLI arguments to the generic scheme.

        This is transitional code while the legacy scheme is still supported.
        """
        translated_args = {"exclude": [], "include": []}

        if self.args.min_date or self.args.max_date:

            translated_args["exclude"].append(
                f"date:"
                + f"min={self.args.min_date}" if self.args.min_date else ""
                + f"max={self.args.max_date}" if self.args.max_date else ""
            )
        if self.args.min_length:
            translated_args["exclude"].append(f"length:min={self.args.min_length}")
        if self.args.non_nucleotide:
            translated_args["exclude"].append("non-nucleotide")
        if self.args.exclude:
            translated_args["exclude"].append(f"name:file={self.args.exclude}")
        if self.args.exclude_where:
            conditions = (
                [self.args.exclude_where]
                if not isinstance(self.args.exclude_where, list)
                else self.args.exclude_where
            )
            for condition in conditions:
                translated_args["exclude"].append(f"metadata:{condition}")
        if self.args.include:
            translated_args["include"].append(f"name:file={self.args.include}")
        if self.args.include_where:
            conditions = (
                [self.args.include_where]
                if not isinstance(self.args.include_where, list)
                else self.args.include_where
            )
            conditions_str = ",".join(conditions)
            translated_args["include"].append(f"metadata:{conditions_str}")

        return translated_args

    @functools.lru_cache()
    def enabled_matchers(self, operation):
        """
        Determine which matchers have been enabled for the given operation.

        Parameters
        ---------
        operation : string
            `"include"` or `"exclude"`

        Returns
        -------
        list
            list of matcher objects that are enabled for the specified operation
        """
        return [
            MATCHER_CLASSES[matcher_type].build(matcher_args)
            for matcher_type, _, matcher_args in (
                filter_directive.partition(":")
                for filter_directive in self.translate_args()[operation]
            )
        ]

    def filter_sequences(self):
        """
        Apply the enabled filters (if any) and apply subsampling (if specified).

        The order of operations is always as follows:
            1. Remove any samples matched by the exclusion filters
            2. Perform subsampling
            3. Add any samples matched by the inclusion filters

        Note that any sample matched by an inclusion rule will _always_ be present in the output.

        Side Effects
        ------------
        self.resulting_sequences : list
            Populated with the result of the specified filtering operation(s)
        """
        self.resulting_sequences = [
            sample for sample in self.sequences if self.no_exclusion_matches(sample)
        ]

        self.resulting_sequences = self.subsample()

        self.resulting_sequences += [
            sample for sample in self.sequences if self.any_inclusion_matches(sample)
        ]

    def no_exclusion_matches(self, sample):
        """
        Determine whether _none_ of the enabled exclusion rules match the given sample.

        Parameters
        ----------
        sample : Sequencse (TODO)
            sample to be evaluated

        Returns
        -------
        bool
            True if none of the enabled exclusion rules match the sample
        """
        return not all(
            [
                filter_obj.is_affected(sample)
                for filter_obj in self.enabled_matchers("exclude")
            ]
        )

    def any_inclusion_matches(self, sample):
        """
        Determine whether _any_ of the specified inclusion rules match the given sample.

        Parameters
        ----------
        sample : Sequencse (TODO)
            sample to be evaluated

        Returns
        -------
        bool
            True if any of the enabled inclusion rules match the sample
        """
        return any(
            [
                filter_obj.is_affected(sample)
                for filter_obj in self.enabled_matchers("include")
            ]
        )

    def subsample(self):
        # TODO
        # Subsampler(self.resulting_sequences, ...).subsample()
        return self.resulting_sequences

    def write_output_file(self):
        num_sequences = len(self.resulting_sequences)
        if num_sequences == 0:
            raise Exception(
                "ERROR: All samples have been dropped! Check filter rules and metadata file format."
            )

        if self.sequence_file.is_vcf:
            VCFSelector(
                selected_samples=self.resulting_sequences,
                samplefile=self.sequence_file,
                output_fname=self.args.output,
            ).write_output_vcf()
        elif self.sequence_file.is_fasta:
            SeqIO.write(
                [sequence.bio_seq for sequence in self.resulting_sequences],
                self.args.output,
                "fasta",
            )
        else:
            raise
        print(f"{num_sequences} sequences written to {self.args.output}")

    def print_summary(self):
        """
        # TODO each filter_obj keeps a count of records affected. or even the actual records, when debug flag.
        total_affected = {"include": 0, "exclude": 0}
        for filter_obj in self.filter_chain:
            op_sign = {"include": "+", "exclude": "-"}[filter_obj.filter_operation()]
            num_affected = len(filter_obj.affected_sequences_set)
            total_affected[filter_obj.filter_operation()] += num_affected
            print(f"{op_sign}{num_affected} samples by {filter_obj.__class__}")

        print(
            "\n"
            f"{total_affected['exclude']} records filtered out."
            f"{total_affected['include']} records re-added."
            f"{total_affected['include'] - total_affected['exclude']} net change."
        )
        """

    def check_vcftools(self):
        if not self.sequence_file.is_vcf:
            # vcftools is irrelevant if the sequences file is not in VCF format
            return

        from shutil import which

        if which("vcftools") is None:
            raise MissingDependencyException("vcftools")
