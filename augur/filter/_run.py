from argparse import Namespace
from tempfile import NamedTemporaryFile
from augur.errors import AugurError

from augur.io.metadata import Metadata
from augur.io.tabular_file import InvalidDelimiter
from . import constants
from .dates import parse_dates
from .debug import print_debug
from .io import get_useful_metadata_columns, initialize_input_source_table, import_metadata, import_sequence_index, print_db_report, write_outputs
from .include_exclude_rules import apply_filters, construct_filters
from .report import print_report
from .subsample import apply_subsampling


def run(args: Namespace):
    with NamedTemporaryFile() as file:
        # Set the database file as a variable that can be easily accessed within
        # functions deep in the call stack. It could be passed down by value,
        # but that would be tedious and makes it harder to trace references back
        # to the source.
        constants.RUNTIME_DB_FILE = file.name
        constants.RUNTIME_DEBUG = args.debug

        print_debug(f"Temporary database file: {constants.RUNTIME_DB_FILE!r}")

        initialize_input_source_table()

        try:
            metadata = Metadata(args.metadata, id_columns=args.metadata_id_columns, delimiters=args.metadata_delimiters)
        except InvalidDelimiter:
            raise AugurError(
                f"Could not determine the delimiter of {args.metadata!r}. "
                f"Valid delimiters are: {args.metadata_delimiters!r}. "
                "This can be changed with --metadata-delimiters."
            )
        columns = get_useful_metadata_columns(args, metadata.id_column, metadata.columns)
        import_metadata(metadata, columns)

        import_sequence_index(args)

        parse_dates()

        exclude_by, include_by = construct_filters(args)
        apply_filters(exclude_by, include_by)

        if args.group_by or args.subsample_max_sequences:
            apply_subsampling(args)

        write_outputs(args)

        print_report(args)

        print_db_report()

    # TODO: The current implementation assumes the database file is hidden from
    # the user. If this ever changes, clean the database of any
    # tables/indexes/etc.
