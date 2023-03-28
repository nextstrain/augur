import argparse
import csv
from argparse import Namespace
import os
import re
from textwrap import dedent
from typing import Sequence, Set
import numpy as np
import pandas as pd
from tempfile import NamedTemporaryFile
from collections import defaultdict
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
from augur.io.vcf import is_vcf, write_vcf
from . import constants
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


def read_priority_scores(fname):
    def constant_factory(value):
        return lambda: value

    try:
        with open_file(fname) as pfile:
            return defaultdict(constant_factory(-np.inf), {
                elems[0]: float(elems[1])
                for elems in (line.strip().split('\t') if '\t' in line else line.strip().split() for line in pfile.readlines())
            })
    except Exception:
        raise AugurError(f"missing or malformed priority scores file {fname}")


def write_metadata_based_outputs(input_metadata_path: str, delimiters: Sequence[str],
                                 id_columns: Sequence[str], output_metadata_path: str,
                                 output_strains_path: str, ids_to_write: Set[str]):
    """
    Write output metadata and/or strains file given input metadata information
    and a set of IDs to write.
    """
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
        constants.sequence_index = pd.read_csv(
            sequence_index_path,
            sep=SEQUENCE_INDEX_DELIMITER,
            index_col=SEQUENCE_INDEX_ID_COLUMN,
            dtype={SEQUENCE_INDEX_ID_COLUMN: "string"},
            **PANDAS_READ_CSV_OPTIONS,
        )

        # Remove temporary index file, if it exists.
        if build_sequence_index:
            os.unlink(sequence_index_path)

        constants.sequence_strains = set(constants.sequence_index.index.values)


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


def read_and_output_sequences(args):
    """Read sequences and output all that passed filtering.
    """
    # Force inclusion of specific strains after filtering and subsampling.
    constants.valid_strains = constants.valid_strains | constants.all_sequences_to_include

    # Write output starting with sequences, if they've been requested. It is
    # possible for the input sequences and sequence index to be out of sync
    # (e.g., the index is a superset of the given sequences input), so we need
    # to update the set of strains to keep based on which strains are actually
    # available.
    if is_vcf(args.sequences):
        if args.output:
            # Get the samples to be deleted, not to keep, for VCF
            dropped_samps = list(constants.sequence_strains - constants.valid_strains)
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

                    if sequence.id in constants.valid_strains:
                        write_sequences(sequence, output_handle, 'fasta')
        else:
            observed_sequence_strains = {sequence.id for sequence in sequences}

        if constants.sequence_strains != observed_sequence_strains:
            # Warn the user if the expected strains from the sequence index are
            # not a superset of the observed strains.
            if constants.sequence_strains is not None and observed_sequence_strains > constants.sequence_strains:
                print_err(
                    "WARNING: The sequence index is out of sync with the provided sequences.",
                    "Metadata and strain output may not match sequence output."
                )

            # Update the set of available sequence strains.
            constants.sequence_strains = observed_sequence_strains


def cleanup_outputs(args):
    """Remove output files. Useful when terminating midway through a loop of metadata chunks."""
    if args.output:
        _try_remove(args.output)
    if args.output_metadata:
        _try_remove(args.output_metadata)
    if args.output_strains:
        _try_remove(args.output_strains)
    if args.output_log:
        _try_remove(args.output_log)


def _try_remove(filepath):
    """Remove a file if it exists."""
    try:
        os.remove(filepath)
    except FileNotFoundError:
        pass
