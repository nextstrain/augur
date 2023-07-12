from argparse import Namespace
import csv
import os
from typing import Set
import numpy as np
from collections import defaultdict

from augur.errors import AugurError
from augur.filter.constants import GROUP_BY_GENERATED_COLUMNS
from augur.io.metadata import METADATA_DATE_COLUMN, Metadata
from .include_exclude_rules import parse_filter_query


def read_priority_scores(fname):
    def constant_factory(value):
        return lambda: value

    try:
        with open(fname, encoding='utf-8') as pfile:
            return defaultdict(constant_factory(-np.inf), {
                elems[0]: float(elems[1])
                for elems in (line.strip().split('\t') if '\t' in line else line.strip().split() for line in pfile.readlines())
            })
    except Exception:
        raise AugurError(f"missing or malformed priority scores file {fname}")


def get_desired_metadata_columns(args: Namespace):
    """Return the column names to be used in augur filter.
    This allows subsetting of the columns on I/O to read only the necessary columns.
    """

    # Build a set of columns.
    columns: Set[str] = set()

    # Get metadata information.
    metadata = Metadata(args.metadata,
        delimiters=args.metadata_delimiters,
        id_columns=args.metadata_id_columns,
    )

    columns.add(metadata.id_column)

    # The date column is always used for date parsing.
    if METADATA_DATE_COLUMN in metadata.columns:
        columns.add(METADATA_DATE_COLUMN)

    if args.group_by:
        # Add columns used for grouping.
        for column in args.group_by:
            if column not in GROUP_BY_GENERATED_COLUMNS:
                columns.add(column)

    if args.exclude_where:
        for exclude_where in args.exclude_where:
            column, op, value = parse_filter_query(exclude_where)
            columns.add(column)

    if args.include_where:
        for include_where in args.include_where:
            column, op, value = parse_filter_query(include_where)
            columns.add(column)

    if args.query:
        for column in extract_variables(args.query):
            columns.add(column)

    return columns


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


import pandas as pd
import re


# from https://stackoverflow.com/a/76536356
def extract_variables(query: str):
    # Track variables in a dictionary to be used as a dictionary of globals.
    variables = {}

    while True:
        try:
            # Try creating a Expr object with the query string and dictionary of globals.
            # This will raise an error as long as the dictionary of globals is incomplete.
            env = pd.core.computation.scope.ensure_scope(level=0, global_dict=variables)
            pd.core.computation.eval.Expr(query, env=env)

            # Exit the loop when evaluation is successful.
            break
        except pd.errors.UndefinedVariableError as e:
            # This relies on the format defined here: https://github.com/pandas-dev/pandas/blob/965ceca9fd796940050d6fc817707bba1c4f9bff/pandas/errors/__init__.py#L401
            name = re.findall("name '(.+?)' is not defined", str(e))[0]

            # Add the name to the globals dictionary with a dummy value.
            variables[name] = None

    return list(variables.keys())


def copy_subset_of_metadata(file_from: str, file_to: str, delimiters, id_columns, strains: Set[str]):
    metadata = Metadata(file_from, delimiters=delimiters, id_columns=id_columns)
    with open(file_to, 'w', newline='') as out:
        writer = csv.writer(out, delimiter=metadata.delimiter, lineterminator=os.linesep)
        writer.writerow(metadata.columns)

        for row in metadata.rows():
            if row[metadata.id_column] in strains:
                writer.writerow([row[column] for column in metadata.columns])
