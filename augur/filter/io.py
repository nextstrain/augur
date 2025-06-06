import argparse
import csv
from argparse import Namespace
import os
import re
from textwrap import dedent
from typing import Sequence, Set
import numpy as np
from collections import defaultdict
from xopen import xopen

from augur.errors import AugurError
from augur.io.file import open_file
from augur.io.metadata import Metadata, METADATA_DATE_COLUMN
from augur.io.print import print_err
from .constants import GROUP_BY_GENERATED_COLUMNS
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
        or (args.group_by and GROUP_BY_GENERATED_COLUMNS.intersection(args.group_by))):
        columns.add(METADATA_DATE_COLUMN)

    if args.group_by:
        group_by_set = set(args.group_by)
        requested_generated_columns = group_by_set & GROUP_BY_GENERATED_COLUMNS

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


def write_output_metadata(input_metadata_path: str, delimiters: Sequence[str],
                          id_columns: Sequence[str], output_metadata_path: str,
                          ids_to_write: Set[str]):
    """
    Write output metadata file given input metadata information and a set of IDs
    to write.
    """
    input_metadata = Metadata(input_metadata_path, delimiters, id_columns)

    with xopen(output_metadata_path, "w", newline="") as output_metadata_handle:
        output_metadata = csv.DictWriter(output_metadata_handle, fieldnames=input_metadata.columns,
                                         delimiter="\t", lineterminator=os.linesep)
        output_metadata.writeheader()

        # Write outputs based on rows in the original metadata.
        for row in input_metadata.rows():
            row_id = row[input_metadata.id_column]
            if row_id in ids_to_write:
                output_metadata.writerow(row)


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


def cleanup_outputs(args):
    """Remove output files. Useful when terminating midway through a loop of metadata chunks."""
    if args.output_sequences:
        _try_remove(args.output_sequences)
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
