import json
import pandas as pd


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
