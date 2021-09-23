"""
Filter and subsample a sequence set.
"""
from Bio import SeqIO
from collections import defaultdict
import csv
import datetime
import heapq
import itertools
import json
import numpy as np
import operator
import os
import pandas as pd
import random
import re
import sys
from tempfile import NamedTemporaryFile
import treetime.utils
from typing import Collection

from .index import index_sequences, index_vcf
from .io import open_file, read_metadata, read_sequences, write_sequences
from .utils import is_vcf as filename_is_vcf, read_vcf, read_strains, get_numerical_dates, run_shell_command, shquote, is_date_ambiguous

comment_char = '#'

SEQUENCE_ONLY_FILTERS = (
    "min_length",
    "non_nucleotide",
)


class FilterException(Exception):
    """Representation of an error that occurred during filtering.
    """
    pass


def write_vcf(input_filename, output_filename, dropped_samps):
    if _filename_gz(input_filename):
        input_arg = "--gzvcf"
    else:
        input_arg = "--vcf"

    if _filename_gz(output_filename):
        output_pipe = "| gzip -c"
    else:
        output_pipe = ""

    drop_args = ["--remove-indv " + shquote(s) for s in dropped_samps]

    call = ["vcftools"] + drop_args + [input_arg, shquote(input_filename), "--recode --stdout", output_pipe, ">", shquote(output_filename)]

    print("Filtering samples using VCFTools with the call:")
    print(" ".join(call))
    run_shell_command(" ".join(call), raise_errors = True)
    # remove vcftools log file
    try:
        os.remove('out.log')
    except OSError:
        pass

def read_priority_scores(fname):
    try:
        with open(fname, encoding='utf-8') as pfile:
            return defaultdict(float, {
                elems[0]: float(elems[1])
                for elems in (line.strip().split('\t') if '\t' in line else line.strip().split() for line in pfile.readlines())
            })
    except Exception as e:
        print(f"ERROR: missing or malformed priority scores file {fname}", file=sys.stderr)
        raise e

# Define metadata filters.

def filter_by_exclude_all(metadata):
    """Exclude all strains regardless of the given metadata content.

    This is a placeholder function that can be called as part of a generalized
    loop through all possible functions.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name

    Returns
    -------
    set[str]:
        Empty set of strains

    >>> metadata = pd.DataFrame([{"region": "Africa"}, {"region": "Europe"}], index=["strain1", "strain2"])
    >>> filter_by_exclude_all(metadata)
    set()
    """
    return set()


def filter_by_exclude(metadata, exclude_file):
    """Exclude the given set of strains from the given metadata.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    exclude_file : str
        Filename with strain names to exclude from the given metadata

    Returns
    -------
    set[str]:
        Strains that pass the filter

    >>> metadata = pd.DataFrame([{"region": "Africa"}, {"region": "Europe"}], index=["strain1", "strain2"])
    >>> with NamedTemporaryFile(delete=False) as exclude_file:
    ...     characters_written = exclude_file.write(b'strain1')
    >>> filter_by_exclude(metadata, exclude_file.name)
    {'strain2'}
    >>> os.unlink(exclude_file.name)
    """
    excluded_strains = read_strains(exclude_file)
    return set(metadata.index.values) - excluded_strains


def parse_filter_query(query):
    """Parse an augur filter-style query and return the corresponding column,
    operator, and value for the query.

    Parameters
    ----------
    query : str
        augur filter-style query following the pattern of `"property=value"` or `"property!=value"`

    Returns
    -------
    str :
        Name of column to query
    callable :
        Operator function to test equality or non-equality of values
    str :
        Value of column to query

    >>> parse_filter_query("property=value")
    ('property', <built-in function eq>, 'value')
    >>> parse_filter_query("property!=value")
    ('property', <built-in function ne>, 'value')

    """
    column, value = re.split(r'!?=', query)
    op = operator.eq
    if "!=" in query:
        op = operator.ne

    return column, op, value


def filter_by_exclude_where(metadata, exclude_where):
    """Exclude all strains from the given metadata that match the given exclusion query.

    Unlike pandas query syntax, exclusion queries should follow the pattern of
    `"property=value"` or `"property!=value"`. Additionally, this filter treats
    all values like lowercase strings, so we convert all values to strings first
    and then lowercase them before testing the given query.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    exclude_where : str
        Filter query used to exclude strains

    Returns
    -------
    set[str]:
        Strains that pass the filter

    >>> metadata = pd.DataFrame([{"region": "Africa"}, {"region": "Europe"}], index=["strain1", "strain2"])
    >>> filter_by_exclude_where(metadata, "region!=Europe")
    {'strain2'}
    >>> filter_by_exclude_where(metadata, "region=Europe")
    {'strain1'}
    >>> filter_by_exclude_where(metadata, "region=europe")
    {'strain1'}

    If the column referenced in the given query does not exist, skip the filter.

    >>> sorted(filter_by_exclude_where(metadata, "missing_column=value"))
    ['strain1', 'strain2']

    """
    column, op, value = parse_filter_query(exclude_where)
    if column in metadata.columns:
        # Apply a test operator (equality or inequality) to values from the
        # column in the given query. This produces an array of boolean values we
        # can index with.
        excluded = op(
            metadata[column].astype(str).str.lower(),
            value.lower()
        )

        # Negate the boolean index of excluded strains to get the index of
        # strains that passed the filter.
        included = ~excluded
        filtered = set(metadata[included].index.values)
    else:
        # Skip the filter, if the requested column does not exist.
        filtered = set(metadata.index.values)

    return filtered


def filter_by_query(metadata, query):
    """Filter metadata in the given pandas DataFrame with a query string and return
    the strain names that pass the filter.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    query : str
        Query string for the dataframe.

    Returns
    -------
    set[str]:
        Strains that pass the filter

    >>> metadata = pd.DataFrame([{"region": "Africa"}, {"region": "Europe"}], index=["strain1", "strain2"])
    >>> filter_by_query(metadata, "region == 'Africa'")
    {'strain1'}
    >>> filter_by_query(metadata, "region == 'North America'")
    set()

    """
    return set(metadata.query(query).index.values)


def filter_by_ambiguous_date(metadata, date_column="date", ambiguity="any"):
    """Filter metadata in the given pandas DataFrame where values in the given date
    column have a given level of ambiguity.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    date_column : str
        Column in the dataframe with dates.
    ambiguity : str
        Level of date ambiguity to filter metadata by

    Returns
    -------
    set[str]:
        Strains that pass the filter

    >>> metadata = pd.DataFrame([{"region": "Africa", "date": "2020-01-XX"}, {"region": "Europe", "date": "2020-01-02"}], index=["strain1", "strain2"])
    >>> filter_by_ambiguous_date(metadata)
    {'strain2'}
    >>> sorted(filter_by_ambiguous_date(metadata, ambiguity="month"))
    ['strain1', 'strain2']

    If the requested date column does not exist, we quietly skip this filter.

    >>> sorted(filter_by_ambiguous_date(metadata, date_column="missing_column"))
    ['strain1', 'strain2']

    """
    if date_column in metadata.columns:
        date_is_ambiguous = metadata[date_column].apply(
            lambda date: is_date_ambiguous(date, ambiguity)
        )
        filtered = set(metadata[~date_is_ambiguous].index.values)
    else:
        filtered = set(metadata.index.values)

    return filtered


def filter_by_date(metadata, date_column="date", min_date=None, max_date=None):
    """Filter metadata by minimum or maximum date.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    date_column : str
        Column in the dataframe with dates.
    min_date : float
        Minimum date
    max_date : float
        Maximum date

    Returns
    -------
    set[str]:
        Strains that pass the filter

    >>> metadata = pd.DataFrame([{"region": "Africa", "date": "2020-01-01"}, {"region": "Europe", "date": "2020-01-02"}], index=["strain1", "strain2"])
    >>> filter_by_date(metadata, min_date=numeric_date("2020-01-02"))
    {'strain2'}
    >>> filter_by_date(metadata, max_date=numeric_date("2020-01-01"))
    {'strain1'}
    >>> filter_by_date(metadata, min_date=numeric_date("2020-01-03"), max_date=numeric_date("2020-01-10"))
    set()
    >>> sorted(filter_by_date(metadata, min_date=numeric_date("2019-12-30"), max_date=numeric_date("2020-01-10")))
    ['strain1', 'strain2']
    >>> sorted(filter_by_date(metadata))
    ['strain1', 'strain2']

    If the requested date column does not exist, we quietly skip this filter.

    >>> sorted(filter_by_date(metadata, date_column="missing_column", min_date=numeric_date("2020-01-02")))
    ['strain1', 'strain2']

    """
    strains = set(metadata.index.values)

    # Skip this filter if no valid min/max date is given or the date column does
    # not exist.
    if (not min_date and not max_date) or date_column not in metadata.columns:
        return strains

    dates = get_numerical_dates(metadata, date_col=date_column, fmt="%Y-%m-%d")
    filtered = {strain for strain in strains if dates[strain] is not None}

    if min_date:
        filtered = {s for s in filtered if (np.isscalar(dates[s]) or all(dates[s])) and np.max(dates[s]) >= min_date}

    if max_date:
        filtered = {s for s in filtered if (np.isscalar(dates[s]) or all(dates[s])) and np.min(dates[s]) <= max_date}

    return filtered


def filter_by_sequence_index(metadata, sequence_index):
    """Filter metadata by presence of corresponding entries in a given sequence
    index. This filter effectively intersects the strain ids in the metadata and
    sequence index.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    sequence_index : pandas.DataFrame
        Sequence index

    Returns
    -------
    set[str]:
        Strains that pass the filter

    >>> metadata = pd.DataFrame([{"region": "Africa", "date": "2020-01-01"}, {"region": "Europe", "date": "2020-01-02"}], index=["strain1", "strain2"])
    >>> sequence_index = pd.DataFrame([{"strain": "strain1", "ACGT": 28000}]).set_index("strain")
    >>> filter_by_sequence_index(metadata, sequence_index)
    {'strain1'}

    """
    metadata_strains = set(metadata.index.values)
    sequence_index_strains = set(sequence_index.index.values)

    return metadata_strains & sequence_index_strains


def filter_by_sequence_length(metadata, sequence_index, min_length=0):
    """Filter metadata by sequence length from a given sequence index.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    sequence_index : pandas.DataFrame
        Sequence index
    min_length : int
        Minimum number of standard nucleotide characters (A, C, G, or T) in each sequence

    Returns
    -------
    set[str]:
        Strains that pass the filter

    >>> metadata = pd.DataFrame([{"region": "Africa", "date": "2020-01-01"}, {"region": "Europe", "date": "2020-01-02"}], index=["strain1", "strain2"])
    >>> sequence_index = pd.DataFrame([{"strain": "strain1", "A": 7000, "C": 7000, "G": 7000, "T": 7000}, {"strain": "strain2", "A": 6500, "C": 6500, "G": 6500, "T": 6500}]).set_index("strain")
    >>> filter_by_sequence_length(metadata, sequence_index, min_length=27000)
    {'strain1'}

    It is possible for the sequence index to be missing strains present in the metadata.

    >>> sequence_index = pd.DataFrame([{"strain": "strain3", "A": 7000, "C": 7000, "G": 7000, "T": 7000}, {"strain": "strain2", "A": 6500, "C": 6500, "G": 6500, "T": 6500}]).set_index("strain")
    >>> filter_by_sequence_length(metadata, sequence_index, min_length=27000)
    set()

    """
    strains = set(metadata.index.values)
    filtered_sequence_index = sequence_index.loc[
        sequence_index.index.intersection(strains)
    ]
    filtered_sequence_index["ACGT"] = filtered_sequence_index.loc[:, ["A", "C", "G", "T"]].sum(axis=1)

    return set(filtered_sequence_index[filtered_sequence_index["ACGT"] >= min_length].index.values)


def filter_by_non_nucleotide(metadata, sequence_index):
    """Filter metadata for strains with invalid nucleotide content.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    sequence_index : pandas.DataFrame
        Sequence index

    Returns
    -------
    set[str]:
        Strains that pass the filter

    >>> metadata = pd.DataFrame([{"region": "Africa", "date": "2020-01-01"}, {"region": "Europe", "date": "2020-01-02"}], index=["strain1", "strain2"])
    >>> sequence_index = pd.DataFrame([{"strain": "strain1", "invalid_nucleotides": 0}, {"strain": "strain2", "invalid_nucleotides": 1}]).set_index("strain")
    >>> filter_by_non_nucleotide(metadata, sequence_index)
    {'strain1'}

    """
    strains = set(metadata.index.values)
    filtered_sequence_index = sequence_index.loc[
        sequence_index.index.intersection(strains)
    ]
    no_invalid_nucleotides = filtered_sequence_index["invalid_nucleotides"] == 0

    return set(filtered_sequence_index[no_invalid_nucleotides].index.values)


def include(metadata, include_file):
    """Include strains in the given text file from the given metadata.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    include_file : str
        Filename with strain names to include from the given metadata

    Returns
    -------
    set[str]:
        Strains that pass the filter

    >>> metadata = pd.DataFrame([{"region": "Africa"}, {"region": "Europe"}], index=["strain1", "strain2"])
    >>> with NamedTemporaryFile(delete=False) as include_file:
    ...     characters_written = include_file.write(b'strain1')
    >>> include(metadata, include_file.name)
    {'strain1'}
    >>> os.unlink(include_file.name)

    """
    included_strains = read_strains(include_file)
    return set(metadata.index.values) & included_strains


def include_by_include_where(metadata, include_where):
    """Include all strains from the given metadata that match the given query.

    Unlike pandas query syntax, inclusion queries should follow the pattern of
    `"property=value"` or `"property!=value"`. Additionally, this filter treats
    all values like lowercase strings, so we convert all values to strings first
    and then lowercase them before testing the given query.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    include_where : str
        Filter query used to include strains

    Returns
    -------
    set[str]:
        Strains that pass the filter

    >>> metadata = pd.DataFrame([{"region": "Africa"}, {"region": "Europe"}], index=["strain1", "strain2"])
    >>> include_by_include_where(metadata, "region!=Europe")
    {'strain1'}
    >>> include_by_include_where(metadata, "region=Europe")
    {'strain2'}
    >>> include_by_include_where(metadata, "region=europe")
    {'strain2'}

    If the column referenced in the given query does not exist, skip the filter.

    >>> include_by_include_where(metadata, "missing_column=value")
    set()

    """
    column, op, value = parse_filter_query(include_where)

    if column in metadata.columns:
        # Apply a test operator (equality or inequality) to values from the
        # column in the given query. This produces an array of boolean values we
        # can index with.
        included_index = op(
            metadata[column].astype(str).str.lower(),
            value.lower()
        )
        included = set(metadata[included_index].index.values)
    else:
        # Skip the inclusion filter if the requested column does not exist.
        included = set()

    return included


def construct_filters(args, sequence_index):
    """Construct lists of filters and inclusion criteria based on user-provided
    arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments provided by the user.
    sequence_index : pandas.DataFrame
        Sequence index for the provided arguments.

    Returns
    -------
    list :
        A list of 2-element tuples with a callable to use as a filter and a
        dictionary of kwargs to pass to the callable.
    list :
        A list of 2-element tuples with a callable and dictionary of kwargs that
        determines whether to force include strains in the final output.

    """
    exclude_by = []
    include_by = []

    # Force include sequences specified in file(s).
    if args.include:
        # Collect the union of all given strains to include.
        for include_file in args.include:
            include_by.append((
                include,
                {
                    "include_file": include_file,
                }
            ))

    # Add sequences with particular metadata attributes.
    if args.include_where:
        for include_where in args.include_where:
            include_by.append((
                include_by_include_where,
                {
                    "include_where": include_where,
                }
            ))

    # Exclude all strains by default.
    if args.exclude_all:
        exclude_by.append((filter_by_exclude_all, {}))

    # Filter by sequence index.
    if sequence_index is not None:
        exclude_by.append((
            filter_by_sequence_index,
            {
                "sequence_index": sequence_index,
            },
        ))

    # Remove strains explicitly excluded by name.
    if args.exclude:
        for exclude_file in args.exclude:
            exclude_by.append((
                filter_by_exclude,
                {
                    "exclude_file": exclude_file,
                }
            ))

    # Exclude strain my metadata field like 'host=camel'.
    if args.exclude_where:
        for exclude_where in args.exclude_where:
            exclude_by.append((
                filter_by_exclude_where,
                {"exclude_where": exclude_where}
            ))

    # Exclude strains by metadata, using pandas querying.
    if args.query:
        exclude_by.append((
            filter_by_query,
            {"query": args.query}
        ))

    # Filter by ambiguous dates.
    if args.exclude_ambiguous_dates_by:
        exclude_by.append((
            filter_by_ambiguous_date,
            {
                "date_column": "date",
                "ambiguity": args.exclude_ambiguous_dates_by,
            }
        ))

    # Filter by date.
    if args.min_date or args.max_date:
        exclude_by.append((
            filter_by_date,
            {
                "date_column": "date",
                "min_date": args.min_date,
                "max_date": args.max_date,
            }
        ))

    # Filter by sequence length.
    if args.min_length:
        # Skip VCF files and warn the user that the min length filter does not
        # make sense for VCFs.
        is_vcf = filename_is_vcf(args.sequences)

        if is_vcf: #doesn't make sense for VCF, ignore.
            print("WARNING: Cannot use min_length for VCF files. Ignoring...")
        else:
            exclude_by.append((
                filter_by_sequence_length,
                {
                    "sequence_index": sequence_index,
                    "min_length": args.min_length,
                }
            ))

    # Exclude sequences with non-nucleotide characters.
    if args.non_nucleotide:
        exclude_by.append((
            filter_by_non_nucleotide,
            {
                "sequence_index": sequence_index,
            }
        ))

    return exclude_by, include_by


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


def apply_filters(metadata, exclude_by, include_by):
    """Apply a list of filters to exclude or force-include records from the given
    metadata and return the strains to keep, to exclude, and to force include.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata to filter
    exclude_by : list[tuple]
        A list of 2-element tuples with a callable to filter by in the first
        index and a dictionary of kwargs to pass to the function in the second
        index.
    include_by : list[tuple]
        A list of 2-element tuples in the same format as the ``exclude_by``
        argument.

    Returns
    -------
    set :
        Strains to keep (those that passed all filters)
    list[dict] :
        Strains to exclude along with the function that filtered them and the arguments used to run the function.
    list[dict] :
        Strains to force-include along with the function that filtered them and the arguments used to run the function.

    For example, filter data by minimum date, but force the include of strains
    from Africa.

    >>> metadata = pd.DataFrame([{"region": "Africa", "date": "2020-01-01"}, {"region": "Europe", "date": "2020-10-02"}, {"region": "North America", "date": "2020-01-01"}], index=["strain1", "strain2", "strain3"])
    >>> exclude_by = [(filter_by_date, {"min_date": numeric_date("2020-04-01")})]
    >>> include_by = [(include_by_include_where, {"include_where": "region=Africa"})]
    >>> strains_to_keep, strains_to_exclude, strains_to_include = apply_filters(metadata, exclude_by, include_by)
    >>> strains_to_keep
    {'strain2'}
    >>> sorted(strains_to_exclude, key=lambda record: record["strain"])
    [{'strain': 'strain1', 'filter': 'filter_by_date', 'kwargs': '[["min_date", 2020.25]]'}, {'strain': 'strain3', 'filter': 'filter_by_date', 'kwargs': '[["min_date", 2020.25]]'}]
    >>> strains_to_include
    [{'strain': 'strain1', 'filter': 'include_by_include_where', 'kwargs': '[["include_where", "region=Africa"]]'}]

    We also want to filter by characteristics of the sequence data that we've
    annotated in a sequence index.

    >>> sequence_index = pd.DataFrame([{"strain": "strain1", "A": 7000, "C": 7000, "G": 7000, "T": 7000}, {"strain": "strain2", "A": 6500, "C": 6500, "G": 6500, "T": 6500}, {"strain": "strain3", "A": 1250, "C": 1250, "G": 1250, "T": 1250}]).set_index("strain")
    >>> exclude_by = [(filter_by_sequence_length, {"sequence_index": sequence_index, "min_length": 27000})]
    >>> include_by = [(include_by_include_where, {"include_where": "region=Europe"})]
    >>> strains_to_keep, strains_to_exclude, strains_to_include = apply_filters(metadata, exclude_by, include_by)
    >>> strains_to_keep
    {'strain1'}
    >>> sorted(strains_to_exclude, key=lambda record: record["strain"])
    [{'strain': 'strain2', 'filter': 'filter_by_sequence_length', 'kwargs': '[["min_length", 27000]]'}, {'strain': 'strain3', 'filter': 'filter_by_sequence_length', 'kwargs': '[["min_length", 27000]]'}]
    >>> strains_to_include
    [{'strain': 'strain2', 'filter': 'include_by_include_where', 'kwargs': '[["include_where", "region=Europe"]]'}]

    """
    strains_to_keep = set(metadata.index.values)
    strains_to_filter = []
    strains_to_force_include = []
    distinct_strains_to_force_include = set()

    # Track strains that should be included regardless of filters.
    for include_function, include_kwargs in include_by:
        passed = metadata.pipe(
            include_function,
            **include_kwargs,
        )
        distinct_strains_to_force_include = distinct_strains_to_force_include | passed

        # Track the reason why strains were included.
        if len(passed) > 0:
            include_name = include_function.__name__
            include_kwargs_str = filter_kwargs_to_str(include_kwargs)
            for strain in passed:
                strains_to_force_include.append({
                    "strain": strain,
                    "filter": include_name,
                    "kwargs": include_kwargs_str,
                })

    for filter_function, filter_kwargs in exclude_by:
        # Apply the current function with its given arguments. Each function
        # returns a set of strains that passed the corresponding filter.
        passed = metadata.pipe(
            filter_function,
            **filter_kwargs,
        )

        # Track the strains that failed this filter, so we can explain why later
        # on and update the list of strains to keep to intersect with the
        # strains that passed.
        failed = strains_to_keep - passed
        strains_to_keep = (strains_to_keep & passed)

        # Track the reason each strain was filtered for downstream reporting.
        if len(failed) > 0:
            # Use a human-readable name for each filter when reporting why a strain
            # was excluded.
            filter_name = filter_function.__name__
            filter_kwargs_str = filter_kwargs_to_str(filter_kwargs)
            for strain in failed:
                strains_to_filter.append({
                    "strain": strain,
                    "filter": filter_name,
                    "kwargs": filter_kwargs_str,
                })

        # Stop applying filters if no strains remain.
        if len(strains_to_keep) == 0:
            break

    return strains_to_keep, strains_to_filter, strains_to_force_include


def get_groups_for_subsampling(strains, metadata, group_by=None):
    """Return a list of groups for each given strain based on the corresponding
    metadata and group by column.

    Parameters
    ----------
    strains : list
        A list of strains to get groups for.
    metadata : pandas.DataFrame
        Metadata to inspect for the given strains.
    group_by : list
        A list of metadata (or calculated) columns to group records by.

    Returns
    -------
    dict :
        A mapping of strain names to tuples corresponding to the values of the strain's group.
    list :
        A list of dictionaries with strains that were skipped from grouping and the reason why (see also: `apply_filters` output).

    >>> strains = ["strain1", "strain2"]
    >>> metadata = pd.DataFrame([{"strain": "strain1", "date": "2020-01-01", "region": "Africa"}, {"strain": "strain2", "date": "2020-02-01", "region": "Europe"}]).set_index("strain")
    >>> group_by = ["region"]
    >>> group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata, group_by)
    >>> group_by_strain
    {'strain1': ('Africa',), 'strain2': ('Europe',)}
    >>> skipped_strains
    []

    If we group by year or month, these groups are calculated from the date
    string.

    >>> group_by = ["year", "month"]
    >>> group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata, group_by)
    >>> group_by_strain
    {'strain1': (2020, (2020, 1)), 'strain2': (2020, (2020, 2))}

    If we omit the grouping columns, the result will group by a dummy column.

    >>> group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata)
    >>> group_by_strain
    {'strain1': ('_dummy',), 'strain2': ('_dummy',)}

    If we try to group by columns that don't exist, we get an error.

    >>> group_by = ["missing_column"]
    >>> get_groups_for_subsampling(strains, metadata, group_by)
    Traceback (most recent call last):
      ...
    augur.filter.FilterException: The specified group-by categories (['missing_column']) were not found. No sequences-per-group sampling will be done.

    If we try to group by some columns that exist and some that don't, we allow
    grouping to continue and print a warning message to stderr.

    >>> group_by = ["year", "month", "missing_column"]
    >>> group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata, group_by)
    >>> group_by_strain
    {'strain1': (2020, (2020, 1), 'unknown'), 'strain2': (2020, (2020, 2), 'unknown')}

    If we group by year month and some records don't have that information in
    their date fields, we should skip those records from the group output and
    track which records were skipped for which reasons.

    >>> metadata = pd.DataFrame([{"strain": "strain1", "date": "", "region": "Africa"}, {"strain": "strain2", "date": "2020-02-01", "region": "Europe"}]).set_index("strain")
    >>> group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata, ["year"])
    >>> group_by_strain
    {'strain2': (2020,)}
    >>> skipped_strains
    [{'strain': 'strain1', 'filter': 'skip_group_by_with_ambiguous_year', 'kwargs': ''}]

    Similarly, if we group by month, we should skip records that don't have
    month information in their date fields.

    >>> metadata = pd.DataFrame([{"strain": "strain1", "date": "2020", "region": "Africa"}, {"strain": "strain2", "date": "2020-02-01", "region": "Europe"}]).set_index("strain")
    >>> group_by_strain, skipped_strains = get_groups_for_subsampling(strains, metadata, ["month"])
    >>> group_by_strain
    {'strain2': ((2020, 2),)}
    >>> skipped_strains
    [{'strain': 'strain1', 'filter': 'skip_group_by_with_ambiguous_month', 'kwargs': ''}]

    """
    if group_by:
        groups = group_by
    else:
        groups = ("_dummy",)

    group_by_strain = {}
    skipped_strains = []
    for strain in strains:
        skip_strain = False
        group = []
        m = metadata.loc[strain].to_dict()
        # collect group specifiers
        for c in groups:
            if c == "_dummy":
                group.append(c)
            elif c in m:
                group.append(m[c])
            elif c in ['month', 'year'] and 'date' in m:
                try:
                    year = int(m["date"].split('-')[0])
                except:
                    skipped_strains.append({
                        "strain": strain,
                        "filter": "skip_group_by_with_ambiguous_year",
                        "kwargs": "",
                    })
                    skip_strain = True
                    break
                if c=='month':
                    try:
                        month = int(m["date"].split('-')[1])
                    except:
                        skipped_strains.append({
                            "strain": strain,
                            "filter": "skip_group_by_with_ambiguous_month",
                            "kwargs": "",
                        })
                        skip_strain = True
                        break

                    group.append((year, month))
                else:
                    group.append(year)
            else:
                group.append('unknown')

        if not skip_strain:
            group_by_strain[strain] = tuple(group)

    # If we could not find any requested categories, we cannot complete subsampling.
    distinct_groups = set(group_by_strain.values())
    if len(distinct_groups) == 1 and ('unknown' in distinct_groups or ('unknown',) in distinct_groups):
        error_message = f"The specified group-by categories ({groups}) were not found. No sequences-per-group sampling will be done."

        if any(x in groups for x in ('year', 'month')):
            error_message += " Note that using 'year' or 'year month' requires a column called 'date'."

        # Raise an exception, since we cannot find the requested groups.
        raise FilterException(error_message)

    # Check to see if some categories are missing to warn the user
    group_by = {
        'date' if cat in ('year', 'month') else cat
        for cat in groups
    }
    missing_cats = [cat for cat in group_by if cat not in metadata.columns.values and cat != "_dummy"]
    if missing_cats:
        error_message = []

        if any(cat != 'date' for cat in missing_cats):
            error_message.append(
                "Some of the specified group-by categories couldn't be found: %s" % ", ".join([str(cat) for cat in missing_cats if cat != 'date'])
            )

        if any(cat == 'date' for cat in missing_cats):
            error_message.append("A 'date' column could not be found to group-by year or month.")

        error_message.append("Filtering by group may behave differently than expected!")

        # Print a warning message, but allow grouping to continue.
        print(
            "WARNING: %s" % "\n".join(error_message),
            file=sys.stderr,
        )

    return group_by_strain, skipped_strains


class PriorityQueue:
    """A priority queue implementation that automatically replaces lower priority
    items in the heap with incoming higher priority items.

    Add a single record to a heap with a maximum of 2 records.

    >>> queue = PriorityQueue(max_size=2)
    >>> queue.add({"strain": "strain1"}, 0.5)
    1

    Add another record with a higher priority. The queue should be at its maximum
    size.

    >>> queue.add({"strain": "strain2"}, 1.0)
    2
    >>> queue.heap
    [(0.5, 0, {'strain': 'strain1'}), (1.0, 1, {'strain': 'strain2'})]
    >>> list(queue.get_items())
    [{'strain': 'strain1'}, {'strain': 'strain2'}]

    Add a higher priority record that causes the queue to exceed its maximum
    size. The resulting queue should contain the two highest priority records
    after the lowest priority record is removed.

    >>> queue.add({"strain": "strain3"}, 2.0)
    2
    >>> list(queue.get_items())
    [{'strain': 'strain2'}, {'strain': 'strain3'}]

    Add a record with the same priority as another record, forcing the duplicate
    to be resolved by removing the oldest entry.

    >>> queue.add({"strain": "strain4"}, 1.0)
    2
    >>> list(queue.get_items())
    [{'strain': 'strain4'}, {'strain': 'strain3'}]

    """
    def __init__(self, max_size):
        """Create a fixed size heap (priority queue)

        """
        self.max_size = max_size
        self.heap = []
        self.counter = itertools.count()

    def add(self, item, priority):
        """Add an item to the queue with a given priority.

        If adding the item causes the queue to exceed its maximum size, replace
        the lowest priority item with the given item. The queue stores items
        with an additional heap id value (a count) to resolve ties between items
        with equal priority (favoring the most recently added item).

        """
        heap_id = next(self.counter)

        if len(self.heap) >= self.max_size:
            heapq.heappushpop(self.heap, (priority, heap_id, item))
        else:
            heapq.heappush(self.heap, (priority, heap_id, item))

        return len(self.heap)

    def get_items(self):
        """Return each item in the queue in order.

        Yields
        ------
        Any
            Item stored in the queue.

        """
        for priority, heap_id, item in self.heap:
            yield item


def create_queues_by_group(groups, max_size, max_attempts=100, random_seed=None):
    """Create a dictionary of priority queues per group for the given maximum size.

    When the maximum size is fractional, probabilistically sample the maximum
    size from a Poisson distribution. Make at least the given number of maximum
    attempts to create queues for which the sum of their maximum sizes is
    greater than zero.

    Create queues for two groups with a fixed maximum size.

    >>> groups = ("2015", "2016")
    >>> queues = create_queues_by_group(groups, 2)
    >>> sum(queue.max_size for queue in queues.values())
    4

    Create queues for two groups with a fractional maximum size. Their total max
    size should still be an integer value greater than zero.

    >>> seed = 314159
    >>> queues = create_queues_by_group(groups, 0.1, random_seed=seed)
    >>> int(sum(queue.max_size for queue in queues.values())) > 0
    True

    A subsequent run of this function with the same groups and random seed
    should produce the same queues and queue sizes.

    >>> more_queues = create_queues_by_group(groups, 0.1, random_seed=seed)
    >>> [queue.max_size for queue in queues.values()] == [queue.max_size for queue in more_queues.values()]
    True

    """
    queues_by_group = {}
    total_max_size = 0
    attempts = 0

    if max_size < 1.0:
        random_generator = np.random.default_rng(random_seed)

    # For small fractional maximum sizes, it is possible to randomly select
    # maximum queue sizes that all equal zero. When this happens, filtering
    # fails unexpectedly. We make multiple attempts to create queues with
    # maximum sizes greater than zero for at least one queue.
    while total_max_size == 0 and attempts < max_attempts:
        for group in sorted(groups):
            if max_size < 1.0:
                queue_max_size = random_generator.poisson(max_size)
            else:
                queue_max_size = max_size

            queues_by_group[group] = PriorityQueue(queue_max_size)

        total_max_size = sum(queue.max_size for queue in queues_by_group.values())
        attempts += 1

    return queues_by_group


def register_arguments(parser):
    input_group = parser.add_argument_group("inputs", "metadata and sequences to be filtered")
    input_group.add_argument('--metadata', required=True, metavar="FILE", help="sequence metadata, as CSV or TSV")
    input_group.add_argument('--sequences', '-s', help="sequences in FASTA or VCF format")
    input_group.add_argument('--sequence-index', help="sequence composition report generated by augur index. If not provided, an index will be created on the fly.")
    input_group.add_argument('--metadata-chunk-size', type=int, default=100000, help="maximum number of metadata records to read into memory at a time. Increasing this number can speed up filtering at the cost of more memory used.")
    input_group.add_argument('--metadata-id-columns', default=["strain", "name"], nargs="+", help="names of valid metadata columns containing identifier information like 'strain' or 'name'")

    metadata_filter_group = parser.add_argument_group("metadata filters", "filters to apply to metadata")
    metadata_filter_group.add_argument(
        '--query',
        help="""Filter samples by attribute.
        Uses Pandas Dataframe querying, see https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query for syntax.
        (e.g., --query "country == 'Colombia'" or --query "(country == 'USA' & (division == 'Washington'))")"""
    )
    metadata_filter_group.add_argument('--min-date', type=numeric_date, help="minimal cutoff for date, the cutoff date is inclusive; may be specified as an Augur-style numeric date (with the year as the integer part) or YYYY-MM-DD")
    metadata_filter_group.add_argument('--max-date', type=numeric_date, help="maximal cutoff for date, the cutoff date is inclusive; may be specified as an Augur-style numeric date (with the year as the integer part) or YYYY-MM-DD")
    metadata_filter_group.add_argument('--exclude-ambiguous-dates-by', choices=['any', 'day', 'month', 'year'],
                                help='Exclude ambiguous dates by day (e.g., 2020-09-XX), month (e.g., 2020-XX-XX), year (e.g., 200X-10-01), or any date fields. An ambiguous year makes the corresponding month and day ambiguous, too, even if those fields have unambiguous values (e.g., "201X-10-01"). Similarly, an ambiguous month makes the corresponding day ambiguous (e.g., "2010-XX-01").')
    metadata_filter_group.add_argument('--exclude', type=str, nargs="+", help="file(s) with list of strains to exclude")
    metadata_filter_group.add_argument('--exclude-where', nargs='+',
                                help="Exclude samples matching these conditions. Ex: \"host=rat\" or \"host!=rat\". Multiple values are processed as OR (matching any of those specified will be excluded), not AND")
    metadata_filter_group.add_argument('--exclude-all', action="store_true", help="exclude all strains by default. Use this with the include arguments to select a specific subset of strains.")
    metadata_filter_group.add_argument('--include', type=str, nargs="+", help="file(s) with list of strains to include regardless of priorities or subsampling")
    metadata_filter_group.add_argument('--include-where', nargs='+',
                                help="Include samples with these values. ex: host=rat. Multiple values are processed as OR (having any of those specified will be included), not AND. This rule is applied last and ensures any sequences matching these rules will be included.")

    sequence_filter_group = parser.add_argument_group("sequence filters", "filters to apply to sequence data")
    sequence_filter_group.add_argument('--min-length', type=int, help="minimal length of the sequences")
    sequence_filter_group.add_argument('--non-nucleotide', action='store_true', help="exclude sequences that contain illegal characters")

    subsample_group = parser.add_argument_group("subsampling", "options to subsample filtered data")
    subsample_group.add_argument('--group-by', nargs='+', help="categories with respect to subsample; two virtual fields, \"month\" and \"year\", are supported if they don't already exist as real fields but a \"date\" field does exist")
    subsample_limits_group = subsample_group.add_mutually_exclusive_group()
    subsample_limits_group.add_argument('--sequences-per-group', type=int, help="subsample to no more than this number of sequences per category")
    subsample_limits_group.add_argument('--subsample-max-sequences', type=int, help="subsample to no more than this number of sequences; can be used without the group_by argument")
    probabilistic_sampling_group = subsample_group.add_mutually_exclusive_group()
    probabilistic_sampling_group.add_argument('--probabilistic-sampling', action='store_true', help="Enable probabilistic sampling during subsampling. This is useful when there are more groups than requested sequences. This option only applies when `--subsample-max-sequences` is provided.")
    probabilistic_sampling_group.add_argument('--no-probabilistic-sampling', action='store_false', dest='probabilistic_sampling')
    subsample_group.add_argument('--priority', type=str, help="""tab-delimited file with list of priority scores for strains (e.g., "<strain>\\t<priority>") and no header.
    When scores are provided, Augur converts scores to floating point values, sorts strains within each subsampling group from highest to lowest priority, and selects the top N strains per group where N is the calculated or requested number of strains per group.
    Higher numbers indicate higher priority.
    Since priorities represent relative values between strains, these values can be arbitrary.""")
    subsample_group.add_argument('--subsample-seed', type=int, help="random number generator seed to allow reproducible subsampling (with same input data).")

    output_group = parser.add_argument_group("outputs", "possible representations of filtered data (at least one required)")
    output_group.add_argument('--output', '--output-sequences', '-o', help="filtered sequences in FASTA format")
    output_group.add_argument('--output-metadata', help="metadata for strains that passed filters")
    output_group.add_argument('--output-strains', help="list of strains that passed filters (no header)")
    output_group.add_argument('--output-log', help="tab-delimited file with one row for each filtered strain and the reason it was filtered. Keyword arguments used for a given filter are reported in JSON format in a `kwargs` column.")

    parser.set_defaults(probabilistic_sampling=True)


def validate_arguments(args):
    """Validate arguments and return a boolean representing whether all validation
    rules succeeded.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments from argparse

    Returns
    -------
    bool :
        Validation succeeded.

    """
    # Don't allow sequence output when no sequence input is provided.
    if args.output and not args.sequences:
        print(
            "ERROR: You need to provide sequences to output sequences.",
            file=sys.stderr)
        return False

    # Confirm that at least one output was requested.
    if not any((args.output, args.output_metadata, args.output_strains)):
        print(
            "ERROR: You need to select at least one output.",
            file=sys.stderr)
        return False

    # Don't allow filtering on sequence-based information, if no sequences or
    # sequence index is provided.
    if not args.sequences and not args.sequence_index and any(getattr(args, arg) for arg in SEQUENCE_ONLY_FILTERS):
        print(
            "ERROR: You need to provide a sequence index or sequences to filter on sequence-specific information.",
            file=sys.stderr)
        return False

    # Set flags if VCF
    is_vcf = filename_is_vcf(args.sequences)

    ### Check users has vcftools. If they don't, a one-blank-line file is created which
    #   allows next step to run but error very badly.
    if is_vcf:
        from shutil import which
        if which("vcftools") is None:
            print("ERROR: 'vcftools' is not installed! This is required for VCF data. "
                  "Please see the augur install instructions to install it.",
                  file=sys.stderr)
            return False

    # If user requested grouping, confirm that other required inputs are provided, too.
    if args.group_by and not any((args.sequences_per_group, args.subsample_max_sequences)):
        print(
            "ERROR: You must specify a number of sequences per group or maximum sequences to subsample.",
            file=sys.stderr
        )
        return False

    return True


def run(args):
    '''
    filter and subsample a set of sequences into an analysis set
    '''
    # Validate arguments before attempting any I/O.
    if not validate_arguments(args):
        return 1

    # Determine whether the sequence index exists or whether should be
    # generated. We need to generate an index if the input sequences are in a
    # VCF, if sequence output has been requested (so we can filter strains by
    # sequences that are present), or if any other sequence-based filters have
    # been requested.
    sequence_strains = None
    sequence_index_path = args.sequence_index
    build_sequence_index = False
    is_vcf = filename_is_vcf(args.sequences)

    if sequence_index_path is None and args.sequences and not args.exclude_all:
        build_sequence_index = True

    if build_sequence_index:
        # Generate the sequence index on the fly, for backwards compatibility
        # with older workflows that don't generate the index ahead of time.
        # Create a temporary index using a random filename to avoid collisions
        # between multiple filter commands.
        with NamedTemporaryFile(delete=False) as sequence_index_file:
            sequence_index_path = sequence_index_file.name

        print(
            "Note: You did not provide a sequence index, so Augur will generate one.",
            "You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.",
            file=sys.stderr
        )

        if is_vcf:
            index_vcf(args.sequences, sequence_index_path)
        else:
            index_sequences(args.sequences, sequence_index_path)

    # Load the sequence index, if a path exists.
    sequence_index = None
    if sequence_index_path:
        sequence_index = pd.read_csv(
            sequence_index_path,
            sep="\t",
            index_col="strain",
        )

        # Remove temporary index file, if it exists.
        if build_sequence_index:
            os.unlink(sequence_index_path)

        # Calculate summary statistics needed for filtering.
        sequence_strains = set(sequence_index.index.values)

    #####################################
    #Filtering steps
    #####################################

    # Setup filters.
    exclude_by, include_by = construct_filters(
        args,
        sequence_index,
    )

    # Setup grouping. We handle the following major use cases:
    #
    # 1. group by and sequences per group defined -> use the given values by the
    # user to identify the highest priority records from each group in a single
    # pass through the metadata.
    #
    # 2. group by and maximum sequences defined -> use the first pass through
    # the metadata to count the number of records in each group, calculate the
    # sequences per group that satisfies the requested maximum, and use a second
    # pass through the metadata to select that many sequences per group.
    #
    # 3. group by not defined but maximum sequences defined -> use a "dummy"
    # group such that we select at most the requested maximum number of
    # sequences in a single pass through the metadata.
    #
    # Each case relies on a priority queue to track the highest priority records
    # per group. In the best case, we can track these records in a single pass
    # through the metadata. In the worst case, we don't know how many sequences
    # per group to use, so we need to calculate this number after the first pass
    # and use a second pass to add records to the queue.
    group_by = args.group_by
    sequences_per_group = args.sequences_per_group
    records_per_group = None

    if group_by and args.subsample_max_sequences:
        # In this case, we need two passes through the metadata with the first
        # pass used to count the number of records per group.
        records_per_group = defaultdict(int)
    elif not group_by and args.subsample_max_sequences:
        group_by = ("_dummy",)
        sequences_per_group = args.subsample_max_sequences

    # If we are grouping data, use queues to store the highest priority strains
    # for each group. When no priorities are provided, they will be randomly
    # generated.
    queues_by_group = None
    if group_by:
        # Use user-defined priorities, if possible. Otherwise, setup a
        # corresponding dictionary that returns a random float for each strain.
        if args.priority:
            priorities = read_priority_scores(args.priority)
        else:
            random_generator = np.random.default_rng(args.subsample_seed)
            priorities = defaultdict(random_generator.random)

    # Setup metadata output. We track whether any records have been written to
    # disk yet through the following variables, to control whether we write the
    # metadata's header and open a new file for writing.
    metadata_header = True
    metadata_mode = "w"

    # Setup strain output.
    if args.output_strains:
        output_strains = open(args.output_strains, "w")

    # Setup logging.
    output_log_writer = None
    if args.output_log:
        # Log the names of strains that were filtered or force-included, so we
        # can properly account for each strain (e.g., including those that were
        # initially filtered for one reason and then included again for another
        # reason).
        output_log = open(args.output_log, "w", newline='')
        output_log_header = ("strain", "filter", "kwargs")
        output_log_writer = csv.DictWriter(
            output_log,
            fieldnames=output_log_header,
            delimiter="\t",
            lineterminator="\n",
        )
        output_log_writer.writeheader()

    # Load metadata. Metadata are the source of truth for which sequences we
    # want to keep in filtered output.
    metadata_strains = set()
    valid_strains = set() # TODO: rename this more clearly
    all_sequences_to_include = set()
    filter_counts = defaultdict(int)

    metadata_reader = read_metadata(
        args.metadata,
        id_columns=args.metadata_id_columns,
        chunk_size=args.metadata_chunk_size,
    )
    for metadata in metadata_reader:
        # Maintain list of all strains seen.
        metadata_strains.update(set(metadata.index.values))

        # Filter metadata.
        seq_keep, sequences_to_filter, sequences_to_include = apply_filters(
            metadata,
            exclude_by,
            include_by,
        )
        valid_strains.update(seq_keep)

        # Track distinct strains to include, so we can write their
        # corresponding metadata, strains, or sequences later, as needed.
        distinct_sequences_to_include = {
            record["strain"]
            for record in sequences_to_include
        }
        all_sequences_to_include.update(distinct_sequences_to_include)

        # Track reasons for filtered or force-included strains, so we can
        # report total numbers filtered and included at the end. Optionally,
        # write out these reasons to a log file.
        for filtered_strain in itertools.chain(sequences_to_filter, sequences_to_include):
            filter_counts[(filtered_strain["filter"], filtered_strain["kwargs"])] += 1

            # Log the names of strains that were filtered or force-included,
            # so we can properly account for each strain (e.g., including
            # those that were initially filtered for one reason and then
            # included again for another reason).
            if args.output_log:
                output_log_writer.writerow(filtered_strain)

        if group_by:
            # If grouping, track the highest priority metadata records or
            # count the number of records per group. First, we need to get
            # the groups for the given records.
            try:
                group_by_strain, skipped_strains = get_groups_for_subsampling(
                    seq_keep,
                    metadata,
                    group_by,
                )

                # Track strains skipped during grouping, so users know why those
                # strains were excluded from the analysis.
                for skipped_strain in skipped_strains:
                    filter_counts[(skipped_strain["filter"], skipped_strain["kwargs"])] += 1
                    valid_strains.remove(skipped_strain["strain"])

                    if args.output_log:
                        output_log_writer.writerow(skipped_strain)

                if args.subsample_max_sequences and records_per_group is not None:
                    # Count the number of records per group. We will use this
                    # information to calculate the number of sequences per group
                    # for the given maximum number of requested sequences.
                    for group in group_by_strain.values():
                        records_per_group[group] += 1
                else:
                    # Track the highest priority records, when we already
                    # know the number of sequences allowed per group.
                    if queues_by_group is None:
                        queues_by_group = {}

                    for strain in sorted(group_by_strain.keys()):
                        # During this first pass, we do not know all possible
                        # groups will be, so we need to build each group's queue
                        # as we first encounter the group.
                        group = group_by_strain[strain]
                        if group not in queues_by_group:
                            queues_by_group[group] = PriorityQueue(
                                max_size=sequences_per_group,
                            )

                        queues_by_group[group].add(
                            metadata.loc[strain],
                            priorities[strain],
                        )
            except FilterException as error:
                # When we cannot group by the requested columns, we print a
                # warning to the user and continue without subsampling or
                # grouping. TODO: We should consider treating this case as
                # an actual error and exiting here with a nonzero code.
                group_by = False
                print(
                    f"WARNING: {error}",
                    file=sys.stderr,
                )

        # Always write out strains that are force-included. Additionally, if
        # we are not grouping, write out metadata and strains that passed
        # filters so far.
        strains_to_write = distinct_sequences_to_include
        if not group_by:
            strains_to_write = strains_to_write | seq_keep

        if args.output_metadata:
            # TODO: wrap logic to write metadata into its own function
            metadata.loc[strains_to_write].to_csv(
                args.output_metadata,
                sep="\t",
                header=metadata_header,
                mode=metadata_mode,
            )
            metadata_header = False
            metadata_mode = "a"

        if args.output_strains:
            # TODO: Output strains will no longer be ordered. This is a
            # small breaking change.
            for strain in strains_to_write:
                output_strains.write(f"{strain}\n")

    # In the worst case, we need to calculate sequences per group from the
    # requested maximum number of sequences and the number of sequences per
    # group. Then, we need to make a second pass through the metadata to find
    # the requested number of records.
    if args.subsample_max_sequences and records_per_group is not None:
        # Calculate sequences per group. If there are more groups than maximum
        # sequences requested, sequences per group will be a floating point
        # value and subsampling will be probabilistic.
        try:
            sequences_per_group = calculate_sequences_per_group(
                args.subsample_max_sequences,
                records_per_group.values(),
                args.probabilistic_sampling,
            )
            print(f"Sampling at {sequences_per_group} per group.")
        except TooManyGroupsError as error:
            print(f"ERROR: {error}", file=sys.stderr)
            sys.exit(1)

        if queues_by_group is None:
            # We know all of the possible groups now from the first pass through
            # the metadata, so we can create queues for all groups at once.
            queues_by_group = create_queues_by_group(
                records_per_group.keys(),
                sequences_per_group,
                random_seed=args.subsample_seed,
            )

        # Make a second pass through the metadata, only considering records that
        # have passed filters.
        metadata_reader = read_metadata(
            args.metadata,
            id_columns=args.metadata_id_columns,
            chunk_size=args.metadata_chunk_size,
        )
        for metadata in metadata_reader:
            # Recalculate groups for subsampling as we loop through the
            # metadata a second time. TODO: We could store these in memory
            # during the first pass, but we want to minimize overall memory
            # usage at the moment.
            seq_keep = set(metadata.index.values) & valid_strains
            group_by_strain, skipped_strains = get_groups_for_subsampling(
                seq_keep,
                metadata,
                group_by,
            )

            for strain in sorted(group_by_strain.keys()):
                group = group_by_strain[strain]
                queues_by_group[group].add(
                    metadata.loc[strain],
                    priorities[strain],
                )

    # If we have any records in queues, we have grouped results and need to
    # stream the highest priority records to the requested outputs.
    num_excluded_subsamp = 0
    if queues_by_group:
        # Populate the set of strains to keep from the records in queues.
        subsampled_strains = set()
        for group, queue in queues_by_group.items():
            records = []
            for record in queue.get_items():
                # Each record is a pandas.Series instance. Track the name of the
                # record, so we can output its sequences later.
                subsampled_strains.add(record.name)

                # Construct a data frame of records to simplify metadata output.
                records.append(record)

                if args.output_strains:
                    # TODO: Output strains will no longer be ordered. This is a
                    # small breaking change.
                    output_strains.write(f"{record.name}\n")

            # Write records to metadata output, if requested.
            if args.output_metadata and len(records) > 0:
                records = pd.DataFrame(records)
                records.to_csv(
                    args.output_metadata,
                    sep="\t",
                    header=metadata_header,
                    mode=metadata_mode,
                )
                metadata_header = False
                metadata_mode = "a"

        # Count and optionally log strains that were not included due to
        # subsampling.
        strains_filtered_by_subsampling = valid_strains - subsampled_strains
        num_excluded_subsamp = len(strains_filtered_by_subsampling)
        if output_log_writer:
            for strain in strains_filtered_by_subsampling:
                output_log_writer.writerow({
                    "strain": strain,
                    "filter": "subsampling",
                    "kwargs": "",
                })

        valid_strains = subsampled_strains

    # Force inclusion of specific strains after filtering and subsampling.
    valid_strains = valid_strains | all_sequences_to_include

    # Write output starting with sequences, if they've been requested. It is
    # possible for the input sequences and sequence index to be out of sync
    # (e.g., the index is a superset of the given sequences input), so we need
    # to update the set of strains to keep based on which strains are actually
    # available.
    if is_vcf:
        if args.output:
            # Get the samples to be deleted, not to keep, for VCF
            dropped_samps = list(sequence_strains - valid_strains)
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

        if sequence_strains != observed_sequence_strains:
            # Warn the user if the expected strains from the sequence index are
            # not a superset of the observed strains.
            if sequence_strains is not None and observed_sequence_strains > sequence_strains:
                print(
                    "WARNING: The sequence index is out of sync with the provided sequences.",
                    "Metadata and strain output may not match sequence output.",
                    file=sys.stderr
                )

            # Update the set of available sequence strains.
            sequence_strains = observed_sequence_strains

    # Calculate the number of strains that don't exist in either metadata or
    # sequences.
    num_excluded_by_lack_of_metadata = 0
    if sequence_strains:
        # Update strains to keep based on available sequence data. This prevents
        # writing out strain lists or metadata for strains that have no
        # sequences.
        valid_strains = valid_strains & sequence_strains

        num_excluded_by_lack_of_metadata = len(sequence_strains - metadata_strains)

    if args.output_strains:
        output_strains.close()

    # Calculate the number of strains passed and filtered.
    total_strains_passed = len(valid_strains)
    total_strains_filtered = len(metadata_strains) + num_excluded_by_lack_of_metadata - total_strains_passed

    print(f"{total_strains_filtered} strains were dropped during filtering")

    if num_excluded_by_lack_of_metadata:
        print(f"\t{num_excluded_by_lack_of_metadata} had no metadata")

    report_template_by_filter_name = {
        "filter_by_sequence_index": "{count} had no sequence data",
        "filter_by_exclude_all": "{count} of these were dropped by `--exclude-all`",
        "filter_by_exclude": "{count} of these were dropped because they were in {exclude_file}",
        "filter_by_exclude_where": "{count} of these were dropped because of '{exclude_where}'",
        "filter_by_query": "{count} of these were filtered out by the query: \"{query}\"",
        "filter_by_ambiguous_date": "{count} of these were dropped because of their ambiguous date in {ambiguity}",
        "filter_by_date": "{count} of these were dropped because of their date (or lack of date)",
        "filter_by_sequence_length": "{count} of these were dropped because they were shorter than minimum length of {min_length}bp",
        "filter_by_non_nucleotide": "{count} of these were dropped because they had non-nucleotide characters",
        "skip_group_by_with_ambiguous_year": "{count} were dropped during grouping due to ambiguous year information",
        "skip_group_by_with_ambiguous_month": "{count} were dropped during grouping due to ambiguous month information",
        "include": "{count} strains were added back because they were in {include_file}",
        "include_by_include_where": "{count} sequences were added back because of '{include_where}'",
    }
    for (filter_name, filter_kwargs), count in filter_counts.items():
        if filter_kwargs:
            parameters = dict(json.loads(filter_kwargs))
        else:
            parameters = {}

        parameters["count"] = count
        print("\t" + report_template_by_filter_name[filter_name].format(**parameters))

    if (group_by and args.sequences_per_group) or args.subsample_max_sequences:
        seed_txt = ", using seed {}".format(args.subsample_seed) if args.subsample_seed else ""
        print("\t%i of these were dropped because of subsampling criteria%s" % (num_excluded_subsamp, seed_txt))

    if total_strains_passed == 0:
        print("ERROR: All samples have been dropped! Check filter rules and metadata file format.", file=sys.stderr)
        return 1

    print(f"{total_strains_passed} strains passed all filters")


def _filename_gz(filename):
    return filename.lower().endswith(".gz")


def numeric_date(date):
    """
    Converts the given *date* string to a :py:class:`float`.

    *date* may be given as a number (a float) with year as the integer part, or
    in the YYYY-MM-DD (ISO 8601) syntax.

    >>> numeric_date("2020.42")
    2020.42
    >>> numeric_date("2020-06-04")
    2020.42486...
    """
    try:
        return float(date)
    except ValueError:
        return treetime.utils.numeric_date(datetime.date(*map(int, date.split("-", 2))))


def calculate_sequences_per_group(target_max_value, counts_per_group, probabilistic=True):
    """Calculate the number of sequences per group for a given maximum number of
    sequences to be returned and the number of sequences in each requested
    group. Optionally, allow the result to be probabilistic such that the mean
    result of a Poisson process achieves the calculated sequences per group for
    the given maximum.

    Parameters
    ----------
    target_max_value : int
        Maximum number of sequences to return by subsampling at some calculated
        number of sequences per group for the given counts per group.
    counts_per_group : list[int]
        A list with the number of sequences in each requested group.
    probabilistic : bool
        Whether to allow probabilistic subsampling when the number of groups
        exceeds the requested maximum.

    Raises
    ------
    TooManyGroupsError :
        When there are more groups than sequences per group and probabilistic
        subsampling is not enabled.

    Returns
    -------
    int or float :
        Number of sequences per group.

    """
    try:
        sequences_per_group = _calculate_sequences_per_group(
            target_max_value,
            counts_per_group,
        )
    except TooManyGroupsError as error:
        if probabilistic:
            sequences_per_group = _calculate_fractional_sequences_per_group(
                target_max_value,
                counts_per_group,
            )
        else:
            raise error

    return sequences_per_group


class TooManyGroupsError(ValueError):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return str(self.msg)


def _calculate_total_sequences(
        hypothetical_spg: float, sequence_lengths: Collection[int],
) -> float:
    # calculate how many sequences we'd keep given a hypothetical spg.
    return sum(
        min(hypothetical_spg, sequence_length)
        for sequence_length in sequence_lengths
    )


def _calculate_sequences_per_group(
        target_max_value: int,
        sequence_lengths: Collection[int]
) -> int:
    """This is partially inspired by
    https://github.com/python/cpython/blob/3.8/Lib/bisect.py

    This should return the spg such that we don't exceed the requested
    number of samples.

    Parameters
    ----------
    target_max_value : int
        the total number of sequences allowed across all groups
    sequence_lengths : Collection[int]
        the number of sequences in each group

    Returns
    -------
    int
        maximum number of sequences allowed per group to meet the required maximum total
        sequences allowed

    >>> _calculate_sequences_per_group(4, [4, 2])
    2
    >>> _calculate_sequences_per_group(2, [4, 2])
    1
    >>> _calculate_sequences_per_group(1, [4, 2])
    Traceback (most recent call last):
        ...
    augur.filter.TooManyGroupsError: Asked to provide at most 1 sequences, but there are 2 groups.
    """

    if len(sequence_lengths) > target_max_value:
        # we have more groups than sequences we are allowed, which is an
        # error.

        raise TooManyGroupsError(
            "Asked to provide at most {} sequences, but there are {} "
            "groups.".format(target_max_value, len(sequence_lengths)))

    lo = 1
    hi = target_max_value

    while hi - lo > 2:
        mid = (hi + lo) // 2
        if _calculate_total_sequences(mid, sequence_lengths) <= target_max_value:
            lo = mid
        else:
            hi = mid

    if _calculate_total_sequences(hi, sequence_lengths) <= target_max_value:
        return int(hi)
    else:
        return int(lo)


def _calculate_fractional_sequences_per_group(
        target_max_value: int,
        sequence_lengths: Collection[int]
) -> float:
    """Returns the fractional sequences per group for the given list of group
    sequences such that the total doesn't exceed the requested number of
    samples.

    Parameters
    ----------
    target_max_value : int
        the total number of sequences allowed across all groups
    sequence_lengths : Collection[int]
        the number of sequences in each group

    Returns
    -------
    float
        fractional maximum number of sequences allowed per group to meet the
        required maximum total sequences allowed

    >>> np.around(_calculate_fractional_sequences_per_group(4, [4, 2]), 4)
    1.9375
    >>> np.around(_calculate_fractional_sequences_per_group(2, [4, 2]), 4)
    0.9688

    Unlike the integer-based version of this function, the fractional version
    can accept a maximum number of sequences that exceeds the number of groups.
    In this case, the function returns a fraction that can be used downstream,
    for example with Poisson sampling.

    >>> np.around(_calculate_fractional_sequences_per_group(1, [4, 2]), 4)
    0.4844
    """
    lo = 1e-5
    hi = target_max_value

    while (hi / lo) > 1.1:
        mid = (lo + hi) / 2
        if _calculate_total_sequences(mid, sequence_lengths) <= target_max_value:
            lo = mid
        else:
            hi = mid

    return (lo + hi) / 2
