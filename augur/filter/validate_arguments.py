from shutil import which
from augur.errors import AugurError
from augur.filter.weights_file import get_weighted_columns
from augur.io.sequences import is_vcf as filename_is_vcf, seqkit


SEQUENCE_ONLY_FILTERS = (
    "min_length",
    "max_length",
    "non_nucleotide",
)


def validate_arguments(args):
    """Validate arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments from argparse
    """
    # Don't allow sequence output when no sequence input is provided.
    if args.output_sequences and not args.sequences:
        raise AugurError("You need to provide sequences to output sequences.")

    # Confirm that at least one output was requested.
    if not any((args.output_sequences, args.output_metadata, args.output_strains)):
        raise AugurError("You need to select at least one output.")

    # Don't allow filtering on sequence-based information, if no sequences or
    # sequence index is provided.
    if not args.sequences and not args.sequence_index and any(getattr(args, arg) for arg in SEQUENCE_ONLY_FILTERS):
        raise AugurError("You need to provide a sequence index or sequences to filter on sequence-specific information.")

    # Confirm that seqkit is installed.
    if args.sequences:
        seqkit()

    # Set flags if VCF
    is_vcf = filename_is_vcf(args.sequences)

    # Confirm that vcftools is installed.
    if is_vcf:
        if which("vcftools") is None:
            raise AugurError("'vcftools' is not installed! This is required for VCF data. "
                  "Please see the augur install instructions to install it.")

    # If user requested grouping, confirm that other required inputs are provided, too.
    if args.group_by and not any((args.sequences_per_group, args.subsample_max_sequences)):
        raise AugurError("You must specify a number of sequences per group or maximum sequences to subsample.")
    
    # Weighted columns must be specified explicitly.
    if args.group_by_weights:
        weighted_columns = get_weighted_columns(args.group_by_weights)
        if (not set(weighted_columns) <= set(args.group_by)):
            raise AugurError("Columns in --group-by-weights must be a subset of columns provided in --group-by.")

    # --output-group-by-sizes is only available for --group-by-weights.
    if args.output_group_by_sizes and not args.group_by_weights:
        raise AugurError(
            "--output-group-by-sizes is only available for --group-by-weights. "
            "It may be added to other sampling methods in the future - see <https://github.com/nextstrain/augur/issues/1590>"
        )

    # --group-by-weights cannot be used with --no-probabilistic-sampling.
    if args.group_by_weights and not args.probabilistic_sampling:
        raise AugurError(
            "--group-by-weights cannot be used with --no-probabilistic-sampling."
        )
