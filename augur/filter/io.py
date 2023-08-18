import argparse
import csv
import os
import re
from typing import Sequence, Set
import numpy as np
from collections import defaultdict
from xopen import xopen

from augur.errors import AugurError
from augur.io.metadata import Metadata


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


def write_metadata_based_outputs(input_metadata_path: str, delimiters: Sequence[str],
                                 id_columns: Sequence[str], output_metadata_path: str,
                                 output_strains_path: str, ids_to_write: Set[str]):
    """
    Write output metadata and/or strains file given input metadata information
    and a set of IDs to write.
    """
    input_metadata = Metadata(input_metadata_path, delimiters, id_columns)

    # Handle all outputs with one pass of metadata. This requires using
    # conditionals both outside of and inside the loop through metadata rows.

    # Make these conditionally set variables available at this scope.
    output_metadata_handle = None
    output_metadata = None
    output_strains = None

    # Set up output streams.
    if output_metadata_path:
        output_metadata_handle = xopen(output_metadata_path, "w")
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
ACCEPTED_TYPES = {'int', 'float', 'str'}

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
