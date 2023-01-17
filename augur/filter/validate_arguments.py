from augur.errors import AugurError
from augur.io.vcf import is_vcf as filename_is_vcf


SEQUENCE_ONLY_FILTERS = (
    "min_length",
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
    if args.output and not args.sequences:
        raise AugurError("You need to provide sequences to output sequences.")

    # Confirm that at least one output was requested.
    if not any((args.output, args.output_metadata, args.output_strains)):
        raise AugurError("You need to select at least one output.")

    # Don't allow filtering on sequence-based information, if no sequences or
    # sequence index is provided.
    if not args.sequences and not args.sequence_index and any(getattr(args, arg) for arg in SEQUENCE_ONLY_FILTERS):
        raise AugurError("You need to provide a sequence index or sequences to filter on sequence-specific information.")

    # Set flags if VCF
    is_vcf = filename_is_vcf(args.sequences)

    # Confirm that vcftools is installed.
    if is_vcf:
        from shutil import which
        if which("vcftools") is None:
            raise AugurError("'vcftools' is not installed! This is required for VCF data. "
                  "Please see the augur install instructions to install it.")

    # If user requested grouping, confirm that other required inputs are provided, too.
    if args.group_by and not any((args.sequences_per_group, args.subsample_max_sequences)):
        raise AugurError("You must specify a number of sequences per group or maximum sequences to subsample.")
