"""
A suite of commands to help with data curation.
"""
import sys
from collections import deque
from textwrap import dedent
from typing import Iterable, Set

from augur.argparse_ import add_command_subparsers
from augur.errors import AugurError
from augur.io.json import dump_ndjson, load_ndjson
from augur.io.metadata import InvalidDelimiter, read_table_to_dict, read_metadata_with_sequences, write_records_to_tsv
from augur.io.sequences import write_records_to_fasta
from . import format_dates, normalize_strings, passthru, titlecase, apply_geolocation_rules, apply_record_annotations, abbreviate_authors, parse_genbank_location, transform_strain_name, rename


SUBCOMMAND_ATTRIBUTE = '_curate_subcommand'
SUBCOMMANDS = [
    passthru,
    normalize_strings,
    format_dates,
    titlecase,
    apply_geolocation_rules,
    apply_record_annotations,
    abbreviate_authors,
    parse_genbank_location,
    transform_strain_name,
    rename,
]


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("curate", help=__doc__)

    # Add print_help so we can run it when no subcommands are called
    parser.set_defaults(print_help = parser.print_help)

    # Add subparsers for subcommands
    subparsers = parser.add_subparsers(dest="subcommand", required=False)
    # Using a subcommand attribute so subcommands are not directly
    # run by top level Augur. Process I/O in `curate`` so individual
    # subcommands do not have to worry about it.
    add_command_subparsers(subparsers, SUBCOMMANDS, SUBCOMMAND_ATTRIBUTE)

    return parser


def validate_records(records: Iterable[dict], subcmd_name: str, is_input: bool) -> Iterable[dict]:
    """
    Validate that the provided *records* all have the same fields.
    Uses the keys of the first record to check against all other records.

    Parameters
    ----------
    records: iterable of dict

    subcmd_name: str
        The name of the subcommand whose output is being validated; used in
        error messages displayed to the user.

    is_input: bool
        Whether the provided records come directly from user provided input
    """
    error_message = "Records do not have the same fields! "
    if is_input:
        error_message += "Please check your input data has the same fields."
    else:
        # Hopefully users should not run into this error as it means we are
        # not uniformly adding/removing fields from records
        error_message += dedent(f"""\
            Something unexpected happened during the augur curate {subcmd_name} command.
            To report this, please open a new issue including the original command:
                <https://github.com/nextstrain/augur/issues/new/choose>
            """)

    first_record_keys: Set[str] = set()
    for idx, record in enumerate(records):
        if idx == 0:
            first_record_keys.update(record.keys())
        else:
            if set(record.keys()) != first_record_keys:
                raise AugurError(error_message)
        yield record


def run(args):
    # Print help if no subcommands are used
    if not getattr(args, SUBCOMMAND_ATTRIBUTE, None):
        return args.print_help()

    # Check provided args are valid and required combination of args are provided
    if not args.fasta and (args.seq_id_column or args.seq_field):
        raise AugurError("The --seq-id-column and --seq-field options should only be used when providing a FASTA file.")

    if args.fasta and (not args.seq_id_column or not args.seq_field):
        raise AugurError("The --seq-id-column and --seq-field options are required for a FASTA file input.")

    if not args.output_fasta and (args.output_id_field or args.output_seq_field):
        raise AugurError("The --output-id-field and --output-seq-field options should only be used when requesting a FASTA output.")

    if args.output_fasta and (not args.output_id_field or not args.output_seq_field):
        raise AugurError("The --output-id-field and --output-seq-field options are required for a FASTA output.")

    # Read inputs
    # Special case single hyphen as stdin
    if args.metadata == '-':
        args.metadata = sys.stdin.buffer

    if args.metadata and args.fasta:
        try:
            records = read_metadata_with_sequences(
                args.metadata,
                args.metadata_delimiters,
                args.fasta,
                args.seq_id_column,
                args.seq_field,
                args.unmatched_reporting,
                args.duplicate_reporting)
        except InvalidDelimiter:
            raise AugurError(
                f"Could not determine the delimiter of {args.metadata!r}. "
                f"Valid delimiters are: {args.metadata_delimiters!r}. "
                "This can be changed with --metadata-delimiters."
            )
    elif args.metadata:
        try:
            records = read_table_to_dict(args.metadata, args.metadata_delimiters, args.duplicate_reporting, args.id_column)
        except InvalidDelimiter:
            raise AugurError(
                f"Could not determine the delimiter of {args.metadata!r}. "
                f"Valid delimiters are: {args.metadata_delimiters!r}. "
                "This can be changed with --metadata-delimiters."
            )
    elif not sys.stdin.isatty():
        records = load_ndjson(sys.stdin)
    else:
        raise AugurError(dedent("""\
            No valid inputs were provided.
            NDJSON records can be streamed from stdin or
            input files can be provided via the command line options `--metadata` and `--fasta`.
            See the command's help message for more details."""))

    # Get the name of the subcmd being run
    subcmd_name = args.subcommand

    # Validate records have the same input fields
    validated_input_records = validate_records(records, subcmd_name, True)

    # Run subcommand to get modified records
    modified_records = getattr(args, SUBCOMMAND_ATTRIBUTE).run(args, validated_input_records)

    # Validate modified records have the same output fields
    validated_output_records = validate_records(modified_records, subcmd_name, False)

    # Output modified records
    # First output FASTA, since the write fasta function yields the records again
    # and removes the sequences from the records
    if args.output_fasta:
        validated_output_records = write_records_to_fasta(
            validated_output_records,
            args.output_fasta,
            args.output_id_field,
            args.output_seq_field)

    if args.output_metadata:
        write_records_to_tsv(validated_output_records, args.output_metadata)

    if not (args.output_fasta or args.output_metadata):
        dump_ndjson(validated_output_records)
    else:
        # Exhaust generator to ensure we run through all records
        # when only a FASTA output is requested but not a metadata output
        deque(validated_output_records, maxlen=0)
