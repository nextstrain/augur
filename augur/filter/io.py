import argparse
import csv
from argparse import Namespace
import os
import re
from textwrap import dedent
from typing import Sequence
from tempfile import NamedTemporaryFile
from xopen import xopen

from augur.errors import AugurError
from augur.index import (
    index_sequences,
    index_vcf,
    ID_COLUMN as SEQUENCE_INDEX_ID_COLUMN,
    DELIMITER as SEQUENCE_INDEX_DELIMITER,
)
from augur.io.file import PANDAS_READ_CSV_OPTIONS, open_file
from augur.io.metadata import Metadata, METADATA_DATE_COLUMN
from augur.io.print import print_err
from augur.io.sequences import read_sequences, write_sequences
from augur.io.sqlite3 import DuplicateError, Sqlite3Database, sanitize_identifier
from augur.io.tabular_file import InvalidDelimiter, TabularFile
from augur.io.vcf import is_vcf, write_vcf
from . import constants
from .debug import print_debug
from .include_exclude_rules import extract_variables, parse_filter_query


def get_useful_metadata_columns(args: Namespace, id_column: str, all_columns: Sequence[str]):
    """Return a list of column names that are used in augur filter.
    This allows reading only the necessary columns.
    """

    # Start with just the ID column.
    columns = {id_column}

    # Add the date column if it is used.
    if (args.exclude_ambiguous_dates_by
        or args.min_date
        or args.max_date
        or (args.group_by and constants.GROUP_BY_GENERATED_COLUMNS.intersection(args.group_by))):
        columns.add(METADATA_DATE_COLUMN)

    if args.group_by:
        group_by_set = set(args.group_by)
        requested_generated_columns = group_by_set & constants.GROUP_BY_GENERATED_COLUMNS

        # Add columns used for grouping.
        columns.update(group_by_set - requested_generated_columns)

        # Show warning for ignored columns.
        ignored_columns = requested_generated_columns.intersection(set(all_columns))
        for col in sorted(ignored_columns):
            print_err(f"WARNING: `--group-by {col}` uses a generated {col} value from the {METADATA_DATE_COLUMN!r} column. The custom '{col}' column in the metadata is ignored for grouping purposes.")

    # Add columns used in exclude queries.
    if args.exclude_where:
        for query in args.exclude_where:
            column, op, value = parse_filter_query(query)
            columns.add(column)

    # Add columns used in include queries.
    if args.include_where:
        for query in args.include_where:
            column, op, value = parse_filter_query(query)
            columns.add(column)

    # Add columns used in Pandas queries.
    if args.query:
        if args.query_columns:
            # Use column names explicitly specified by the user.
            for column, dtype in args.query_columns:
                columns.add(column)

        # Attempt to automatically extract columns from the query.
        variables = extract_variables(args.query)
        if variables is None and not args.query_columns:
            print_err(dedent(f"""\
                WARNING: Could not infer columns from the pandas query. Reading all metadata columns,
                which may impact execution time. If the query is valid, please open a new issue:

                    <https://github.com/nextstrain/augur/issues/new/choose>

                and add the query to the description:

                    {args.query}"""))
            columns.update(all_columns)
        else:
            columns.update(variables)

    return list(columns)


def import_priorities_table(path):
    """Import a priorities file into the database."""
    try:
        priorities = TabularFile(path, delimiters=['\t'], header=False,
                                 columns=[constants.ID_COLUMN, constants.PRIORITY_COLUMN])
        with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
            _import_tabular_file(priorities, db, constants.PRIORITIES_TABLE)
    except (FileNotFoundError, InvalidDelimiter):
        raise AugurError(f"missing or malformed priority scores file {path}")

    try:
        _validate_priorities_table()
    except ValueError:
        # TODO: Surface the underlying error message.
        raise AugurError(f"missing or malformed priority scores file {path}")


def _validate_priorities_table():
    """Query the priorities table and error upon any invalid scores."""
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        result = db.connection.execute(f"""
            SELECT {constants.ID_COLUMN}, {constants.PRIORITY_COLUMN}
            FROM {constants.PRIORITIES_TABLE}
        """)
        for row in result:
            try:
                float(row[constants.PRIORITY_COLUMN])
            except ValueError:
                raise ValueError(f"Priority score for strain '{row[constants.ID_COLUMN]}' ('{row[constants.PRIORITY_COLUMN]}') is not a valid number.")


def _write_metadata_based_outputs(input_metadata_path: str, delimiters: Sequence[str],
                                 id_columns: Sequence[str], output_metadata_path: str,
                                 output_strains_path: str):
    """
    Write output metadata and/or strains file given input metadata information
    and a set of IDs to write.
    """
    ids_to_write = _get_valid_strains()

    input_metadata = Metadata(input_metadata_path, id_columns, delimiters=delimiters)

    # Handle all outputs with one pass of metadata. This requires using
    # conditionals both outside of and inside the loop through metadata rows.

    # Make these conditionally set variables available at this scope.
    output_metadata_handle = None
    output_metadata = None
    output_strains = None

    # Set up output streams.
    if output_metadata_path:
        output_metadata_handle = xopen(output_metadata_path, "w", newline="")
        output_metadata = csv.DictWriter(output_metadata_handle, fieldnames=input_metadata.columns,
                                         delimiter="\t", lineterminator=os.linesep)
        output_metadata.writeheader()
    if output_strains_path:
        output_strains = open(output_strains_path, "w")

    # Write outputs based on rows in the original metadata.
    for row in input_metadata.rows():
        row_id = row[input_metadata.id_column]
        if row_id in ids_to_write:
            if output_metadata:
                output_metadata.writerow(row)
            if output_strains:
                output_strains.write(row_id + '\n')

    # Close file handles.
    if output_metadata_handle:
        output_metadata_handle.close()
    if output_strains:
        output_strains.close()


# These are the types accepted in the following function.
ACCEPTED_TYPES = {'int', 'float', 'bool', 'str'}

def column_type_pair(input: str):
    """Get a 2-tuple for column name to type.

    Intended to be used as the argument type converter for argparse options that
    take type maps in a 'column:type' format.
    """

    match = re.match(f"^(.+?):({'|'.join(ACCEPTED_TYPES)})$", input)
    if not match:
        raise argparse.ArgumentTypeError(f"Column data types must be in the format 'column:type', where type is one of ({','.join(ACCEPTED_TYPES)}).")

    column = match[1]
    dtype = match[2]

    return (column, dtype)


def initialize_input_source_table():
    """Create the input source table without any rows."""
    with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
        db.connection.execute(f"""
            CREATE TABLE {constants.INPUT_SOURCE_TABLE} (
                {constants.ID_COLUMN} TEXT,
                {constants.STRAIN_IN_METADATA_COLUMN} INTEGER,
                {constants.STRAIN_IN_SEQUENCE_INDEX_COLUMN} INTEGER,
                {constants.STRAIN_IN_SEQUENCES_COLUMN} INTEGER
            )
        """)
        db.create_primary_index(constants.INPUT_SOURCE_TABLE, constants.ID_COLUMN)


def _import_tabular_file(file: TabularFile, db: Sqlite3Database, table: str, columns=None):
    """Import a tabular file into a new table in an existing database.

    Parameters
    ----------
    file
        File to import from.
    db
        Database to import into.
    table
        Table name to import into.
    columns
        Columns to import.
    """
    if columns is None:
        columns = file.columns
    else:
        for column in list(columns):
            if column not in file.columns:
                # Ignore missing columns. Don't error since augur filter's
                # --exclude-where allows invalid columns to be specified (they
                # are just ignored).
                print_err(f"WARNING: Column '{column}' does not exist in the metadata file. This may cause subsequent errors.")
                columns.remove(column)
    db.create_table(table, columns)
    db.insert(table, columns, file.rows())


def import_metadata(metadata: Metadata, columns):
    """Import metadata into the database."""
    with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
        _import_tabular_file(metadata, db, constants.METADATA_TABLE, columns)

        try:
            db.create_primary_index(constants.METADATA_TABLE, metadata.id_column)
        except DuplicateError as error:
            duplicates = error.duplicated_values
            raise AugurError(f"The following strains are duplicated in '{metadata.path}':\n" + "\n".join(sorted(duplicates)))

        # If the strain is already in the input source table, update it.
        db.connection.execute(f"""
            UPDATE {constants.INPUT_SOURCE_TABLE}
            SET {constants.STRAIN_IN_METADATA_COLUMN} = TRUE
            WHERE {constants.ID_COLUMN} IN (
                SELECT {sanitize_identifier(metadata.id_column)}
                FROM {constants.METADATA_TABLE}
            )
        """)

        # Otherwise, add an entry.
        db.connection.execute(f"""
            INSERT OR IGNORE INTO {constants.INPUT_SOURCE_TABLE} (
                {constants.ID_COLUMN},
                {constants.STRAIN_IN_METADATA_COLUMN},
                {constants.STRAIN_IN_SEQUENCE_INDEX_COLUMN},
                {constants.STRAIN_IN_SEQUENCES_COLUMN}
            )
            SELECT
                {sanitize_identifier(metadata.id_column)} AS {constants.ID_COLUMN},
                TRUE AS {constants.STRAIN_IN_METADATA_COLUMN},
                FALSE AS {constants.STRAIN_IN_SEQUENCE_INDEX_COLUMN},
                FALSE AS {constants.STRAIN_IN_SEQUENCES_COLUMN}
            FROM {constants.METADATA_TABLE}
        """)


def import_sequence_index(args):
    # Determine whether the sequence index exists or whether should be
    # generated. We need to generate an index if the input sequences are in a
    # VCF, if sequence output has been requested (so we can filter strains by
    # sequences that are present), or if any other sequence-based filters have
    # been requested.
    sequence_index_path = args.sequence_index
    build_sequence_index = False

    # Don't build sequence index with --exclude-all since the only way to add
    # strains back in with this flag are the `--include` or `--include-where`
    # options, so we know we don't need a sequence index to apply any additional
    # filters.
    if sequence_index_path is None and args.sequences and not args.exclude_all:
        build_sequence_index = True

    if build_sequence_index:
        sequence_index_path = _generate_sequence_index(args.sequences)

    # Load the sequence index, if a path exists.
    if sequence_index_path:
        try:
            sequence_index = TabularFile(sequence_index_path, header=True, delimiters=[SEQUENCE_INDEX_DELIMITER])
        except InvalidDelimiter:
            # This can happen for single-column files (e.g. VCF sequence indexes).
            # If so, use a tab character as an arbitrary delimiter.
            sequence_index = TabularFile(sequence_index_path, header=True, delimiter='\t')
        with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
            # Import the sequence index.
            _import_tabular_file(sequence_index, db, constants.SEQUENCE_INDEX_TABLE)
            # FIXME: set type affinity of SEQUENCE_INDEX_ID_COLUMN to TEXT
            db.create_primary_index(constants.SEQUENCE_INDEX_TABLE, SEQUENCE_INDEX_ID_COLUMN)

            # If the strain is already in the input source table, update it.
            db.connection.execute(f"""
                UPDATE {constants.INPUT_SOURCE_TABLE}
                SET {constants.STRAIN_IN_SEQUENCE_INDEX_COLUMN} = TRUE
                WHERE {constants.ID_COLUMN} IN (
                    SELECT {SEQUENCE_INDEX_ID_COLUMN}
                    FROM {constants.SEQUENCE_INDEX_TABLE}
                )
            """)

            # Otherwise, add an entry.
            db.connection.execute(f"""
                INSERT OR IGNORE INTO {constants.INPUT_SOURCE_TABLE} (
                    {constants.ID_COLUMN},
                    {constants.STRAIN_IN_METADATA_COLUMN},
                    {constants.STRAIN_IN_SEQUENCE_INDEX_COLUMN},
                    {constants.STRAIN_IN_SEQUENCES_COLUMN}
                )
                SELECT
                    {SEQUENCE_INDEX_ID_COLUMN} AS {constants.ID_COLUMN},
                    FALSE AS {constants.STRAIN_IN_METADATA_COLUMN},
                    TRUE AS {constants.STRAIN_IN_SEQUENCE_INDEX_COLUMN},
                    FALSE AS {constants.STRAIN_IN_SEQUENCES_COLUMN}
                FROM {constants.SEQUENCE_INDEX_TABLE}
            """)

        # Remove temporary index file, if it exists.
        if build_sequence_index:
            os.unlink(sequence_index_path)


def _generate_sequence_index(sequences_file):
    """Generate a sequence index file.
    """
    # Generate the sequence index on the fly, for backwards compatibility
    # with older workflows that don't generate the index ahead of time.
    # Create a temporary index using a random filename to avoid collisions
    # between multiple filter commands.
    with NamedTemporaryFile(delete=False) as sequence_index_file:
        sequence_index_path = sequence_index_file.name

    print_err(
        "Note: You did not provide a sequence index, so Augur will generate one.",
        "You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`."
    )

    # FIXME: call a function in index_sequences which already handles VCF vs. FASTA
    if is_vcf(sequences_file):
        index_vcf(sequences_file, sequence_index_path)
    else:
        index_sequences(sequences_file, sequence_index_path)

    return sequence_index_path


def write_outputs(args):
    """Write the output files that were requested."""

    _read_and_output_sequences(args)

    _write_metadata_based_outputs(args.metadata, args.metadata_delimiters, args.metadata_id_columns, args.output_metadata, args.output_strains)

    if args.output_log:
        _output_log(args.output_log)


def _read_and_output_sequences(args):
    """Read sequences and output all that passed filtering.
    """
    valid_strains = _get_valid_strains()

    # Write output starting with sequences, if they've been requested. It is
    # possible for the input sequences and sequence index to be out of sync
    # (e.g., the index is a superset of the given sequences input), so we need
    # to update the set of strains to keep based on which strains are actually
    # available.
    if is_vcf(args.sequences):
        if args.output:
            # Get the samples to be deleted, not to keep, for VCF
            dropped_samps = _get_strains_to_drop_from_vcf()
            write_vcf(args.sequences, args.output, dropped_samps)
    elif args.sequences:
        sequences = read_sequences(args.sequences)

        # If the user requested sequence output, stream to disk all sequences
        # that passed all filters to avoid reading sequences into memory first.
        # Even if we aren't emitting sequences, we track the observed strain
        # names in the sequence file as part of the single pass to allow
        # comparison with the provided sequence index.
        if args.output:
            observed_sequence_strains = set()
            with open_file(args.output, "wt") as output_handle:
                for sequence in sequences:
                    observed_sequence_strains.add(sequence.id)

                    if sequence.id in valid_strains:
                        write_sequences(sequence, output_handle, 'fasta')
        else:
            observed_sequence_strains = {sequence.id for sequence in sequences}

        # Update the input source table.
        with Sqlite3Database(constants.RUNTIME_DB_FILE, mode="rw") as db:
            # If the strain is already in the input source table, update it.
            quoted_strains = (f"'{strain}'" for strain in observed_sequence_strains)
            db.connection.execute(f"""
                UPDATE {constants.INPUT_SOURCE_TABLE}
                SET {constants.STRAIN_IN_SEQUENCES_COLUMN} = TRUE
                WHERE {constants.ID_COLUMN} IN ({','.join(quoted_strains)})
            """)

            # Otherwise, add an entry.
            rows = ({'strain': strain} for strain in observed_sequence_strains)
            db.connection.executemany(f"""
                INSERT OR IGNORE INTO {constants.INPUT_SOURCE_TABLE} (
                    {constants.ID_COLUMN},
                    {constants.STRAIN_IN_METADATA_COLUMN},
                    {constants.STRAIN_IN_SEQUENCE_INDEX_COLUMN},
                    {constants.STRAIN_IN_SEQUENCES_COLUMN}
                )
                VALUES (
                    :strain,
                    FALSE,
                    FALSE,
                    TRUE
                )
            """, rows)

        with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
            # Only run this if the sequence index table exists.
            if constants.SEQUENCE_INDEX_TABLE in db.tables():
                result = db.connection.execute(f"""
                    SELECT COUNT(*)
                    FROM {constants.INPUT_SOURCE_TABLE}
                    WHERE {constants.STRAIN_IN_SEQUENCES_COLUMN} AND NOT {constants.STRAIN_IN_SEQUENCE_INDEX_COLUMN}
                """)
                sequences_missing_from_index = result.fetchone()[0]

                if sequences_missing_from_index > 0:
                    # Warn the user if the expected strains from the sequence index are
                    # not a superset of the observed strains.
                    print_err(
                        "WARNING: The sequence index is out of sync with the provided sequences.",
                        "Metadata and strain output may not match sequence output."
                    )


def _output_log(path):
    """Write a file explaining the reason for excluded or force-included strains.

    This file has the following columns:
    1. Strain column
    2. Name of the filter function responsible for inclusion/exclusion
    3. Arguments given to the filter function
    """
    query = f"""
        SELECT
            {constants.ID_COLUMN},
            {constants.FILTER_REASON_COLUMN},
            {constants.FILTER_REASON_KWARGS_COLUMN}
        FROM {constants.FILTER_REASON_TABLE}
        WHERE {constants.FILTER_REASON_COLUMN} IS NOT NULL
    """
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        db.query_to_file(
            query=query,
            path=path,
            header=True,
        )


def _get_valid_strains():
    """Returns the strains that pass all filter rules.
    """
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        result = db.connection.execute(f"""
            SELECT {constants.ID_COLUMN}
            FROM {constants.FILTER_REASON_TABLE}
            WHERE NOT {constants.EXCLUDE_COLUMN} OR {constants.INCLUDE_COLUMN}
        """)
        return {str(row[constants.ID_COLUMN]) for row in result}


def _get_strains_to_drop_from_vcf():
    """Return a set of all strain names that are in the sequence index and did
    not pass filtering and subsampling."""
    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        # Query = all sequence index strains - (all strains that passed).
        # This includes strains that are not present in the metadata.
        result = db.connection.execute(f"""
            SELECT {SEQUENCE_INDEX_ID_COLUMN}
            FROM {constants.SEQUENCE_INDEX_TABLE}
            WHERE {SEQUENCE_INDEX_ID_COLUMN} NOT IN (
                SELECT {constants.ID_COLUMN}
                FROM {constants.FILTER_REASON_TABLE}
                WHERE NOT {constants.EXCLUDE_COLUMN} OR {constants.INCLUDE_COLUMN}
            )
        """)
        return {str(row[SEQUENCE_INDEX_ID_COLUMN]) for row in result}


def print_db_report():
    if not constants.RUNTIME_DEBUG:
        return

    with Sqlite3Database(constants.RUNTIME_DB_FILE) as db:
        result = db.connection.execute(f"""
            SELECT
                name,
                SUM(pgsize) AS size
            FROM dbstat
            GROUP BY name;
        """)
        rows = result.fetchall()

    print_debug(f'The total size of the database was {_human_readable_size(sum(int(row["size"]) for row in rows))}. Breakdown:')

    for row in sorted(rows, key=lambda row: int(row["size"]), reverse=True):
        print_debug(f'{_human_readable_size(row["size"]): >10}  {row["name"]}')


def _human_readable_size(bytes: int, decimal_places=1):
    """Return size in bytes as a human-readable string using larger units.

    Adapted from https://stackoverflow.com/a/43690506
    """
    size = float(bytes)
    units = ['B', 'KiB', 'MiB', 'GiB', 'TiB', 'PiB']
    for unit in units:
        if size < 1024.0 or unit == units[-1]:
            break
        size /= 1024.0
    return f"{size:.{decimal_places}f} {unit}"
