import json
import os
import numpy as np
import pandas as pd
from collections import defaultdict

from augur.errors import AugurError


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


def filter_kwargs_to_str(kwargs):
    """Convert a dictionary of kwargs to a JSON string for downstream reporting.

    This structured string can be converted back into a Python data structure
    later for more sophisticated reporting by specific kwargs.

    This function excludes data types from arguments like pandas DataFrames and
    also converts floating point numbers to a fixed precision for better
    readability and reproducibility.

    Parameters
    ----------
    kwargs : dict
        Dictionary of kwargs passed to a given filter function.

    Returns
    -------
    str :
        String representation of the kwargs for reporting.


    >>> from augur.dates import numeric_date
    >>> from augur.filter.include_exclude_rules import filter_by_sequence_length, filter_by_date
    >>> sequence_index = pd.DataFrame([{"strain": "strain1", "ACGT": 28000}, {"strain": "strain2", "ACGT": 26000}, {"strain": "strain3", "ACGT": 5000}]).set_index("strain")
    >>> exclude_by = [(filter_by_sequence_length, {"sequence_index": sequence_index, "min_length": 27000})]
    >>> filter_kwargs_to_str(exclude_by[0][1])
    '[["min_length", 27000]]'
    >>> exclude_by = [(filter_by_date, {"max_date": numeric_date("2020-04-01"), "min_date": numeric_date("2020-03-01")})]
    >>> filter_kwargs_to_str(exclude_by[0][1])
    '[["max_date", 2020.25], ["min_date", 2020.17]]'

    """
    # Sort keys prior to processing to guarantee the same output order
    # regardless of the input order.
    sorted_keys = sorted(kwargs.keys())

    kwarg_list = []
    for key in sorted_keys:
        value = kwargs[key]

        # Handle special cases for data types that we want to represent
        # differently from their defaults or not at all.
        if isinstance(value, pd.DataFrame):
            continue
        elif isinstance(value, float):
            value = round(value, 2)

        kwarg_list.append((key, value))

    return json.dumps(kwarg_list)


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
