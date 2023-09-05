import json
import operator
import re
import numpy as np
import pandas as pd
from typing import Any, Callable, Dict, List, Set, Tuple

from augur.dates import is_date_ambiguous, get_numerical_dates
from augur.errors import AugurError
from augur.io.metadata import METADATA_DATE_COLUMN
from augur.io.print import print_err
from augur.io.vcf import is_vcf as filename_is_vcf
from augur.utils import read_strains
from . import constants

try:
    # pandas ≥1.5.0 only
    PandasUndefinedVariableError = pd.errors.UndefinedVariableError  # type: ignore
except AttributeError:
    PandasUndefinedVariableError = pd.core.computation.ops.UndefinedVariableError  # type: ignore


# The strains to keep as a result of applying a filter function.
FilterFunctionReturn = Set[str]

# A function to use for filtering. Parameters vary.
FilterFunction = Callable[..., FilterFunctionReturn]

# A dictionary of parameters to pass to the filter function.
FilterFunctionKwargs = Dict[str, Any]

# A pair of values to represent the filter function and parameters to use.
FilterOption = Tuple[FilterFunction, FilterFunctionKwargs]


def filter_by_exclude_all(metadata) -> FilterFunctionReturn:
    """Exclude all strains regardless of the given metadata content.

    This is a placeholder function that can be called as part of a generalized
    loop through all possible functions.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name

    Examples
    --------
    >>> metadata = pd.DataFrame([{"region": "Africa"}, {"region": "Europe"}], index=["strain1", "strain2"])
    >>> filter_by_exclude_all(metadata)
    set()
    """
    return set()


def filter_by_exclude(metadata, exclude_file) -> FilterFunctionReturn:
    """Exclude the given set of strains from the given metadata.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    exclude_file : str
        Filename with strain names to exclude from the given metadata

    Examples
    --------
    >>> import os
    >>> from tempfile import NamedTemporaryFile
    >>> metadata = pd.DataFrame([{"region": "Africa"}, {"region": "Europe"}], index=["strain1", "strain2"])
    >>> with NamedTemporaryFile(delete=False) as exclude_file:
    ...     characters_written = exclude_file.write(b'strain1')
    >>> filter_by_exclude(metadata, exclude_file.name)
    {'strain2'}
    >>> os.unlink(exclude_file.name)
    """
    excluded_strains = read_strains(exclude_file)
    return set(metadata.index.values) - excluded_strains


def _parse_filter_query(query):
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

    Examples
    --------
    >>> _parse_filter_query("property=value")
    ('property', <built-in function eq>, 'value')
    >>> _parse_filter_query("property!=value")
    ('property', <built-in function ne>, 'value')

    """
    column, value = re.split(r'!?=', query)
    op = operator.eq
    if "!=" in query:
        op = operator.ne

    return column, op, value


def filter_by_exclude_where(metadata, exclude_where) -> FilterFunctionReturn:
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

    Examples
    --------
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
    column, op, value = _parse_filter_query(exclude_where)
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


def filter_by_query(metadata, query) -> FilterFunctionReturn:
    """Filter metadata in the given pandas DataFrame with a query string and return
    the strain names that pass the filter.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    query : str
        Query string for the dataframe.

    Examples
    --------
    >>> metadata = pd.DataFrame([{"region": "Africa"}, {"region": "Europe"}], index=["strain1", "strain2"])
    >>> filter_by_query(metadata, "region == 'Africa'")
    {'strain1'}
    >>> filter_by_query(metadata, "region == 'North America'")
    set()

    """
    # Create a copy to prevent modification of the original DataFrame.
    metadata_copy = metadata.copy()

    # Support numeric comparisons in query strings.
    #
    # The built-in data type inference when loading the DataFrame does not
    # support nullable numeric columns, so numeric comparisons won't work on
    # those columns. pd.to_numeric does proper conversion on those columns, and
    # will not make any changes to columns with other values.
    #
    # TODO: Parse the query string and apply conversion only to columns used for
    # numeric comparison. Pandas does not expose the API used to parse the query
    # string internally, so this is non-trivial and requires a bit of
    # reverse-engineering. Commit 2ead5b3e3306dc1100b49eb774287496018122d9 got
    # halfway there but had issues so it was reverted.
    #
    # TODO: Try boolean conversion?
    for column in metadata_copy.columns:
        metadata_copy[column] = pd.to_numeric(metadata_copy[column], errors='ignore')

    return set(metadata_copy.query(query).index.values)


def filter_by_ambiguous_date(metadata, date_column, ambiguity) -> FilterFunctionReturn:
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

    Examples
    --------
    >>> metadata = pd.DataFrame([{"region": "Africa", "date": "2020-01-XX"}, {"region": "Europe", "date": "2020-01-02"}], index=["strain1", "strain2"])
    >>> filter_by_ambiguous_date(metadata, date_column="date", ambiguity="any")
    {'strain2'}
    >>> sorted(filter_by_ambiguous_date(metadata, date_column="date", ambiguity="month"))
    ['strain1', 'strain2']

    If the requested date column does not exist, we quietly skip this filter.

    >>> sorted(filter_by_ambiguous_date(metadata, date_column="missing_column", ambiguity="any"))
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


def skip_group_by_with_ambiguous_year(metadata, date_column) -> FilterFunctionReturn:
    """Alias to filter_by_ambiguous_date for year. This is to have a named function available for the filter reason."""
    return filter_by_ambiguous_date(metadata, date_column, ambiguity="year")


def skip_group_by_with_ambiguous_month(metadata, date_column) -> FilterFunctionReturn:
    """Alias to filter_by_ambiguous_date for month. This is to have a named function available for the filter reason."""
    return filter_by_ambiguous_date(metadata, date_column, ambiguity="month")


def skip_group_by_with_ambiguous_day(metadata, date_column) -> FilterFunctionReturn:
    """Alias to filter_by_ambiguous_date for day. This is to have a named function available for the filter reason."""
    return filter_by_ambiguous_date(metadata, date_column, ambiguity="day")


def filter_by_min_date(metadata, date_column, min_date) -> FilterFunctionReturn:
    """Filter metadata by minimum date.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    date_column : str
        Column in the dataframe with dates.
    min_date : float
        Minimum date

    Examples
    --------
    >>> from augur.dates import numeric_date
    >>> metadata = pd.DataFrame([{"region": "Africa", "date": "2020-01-01"}, {"region": "Europe", "date": "2020-01-02"}], index=["strain1", "strain2"])
    >>> filter_by_min_date(metadata, date_column="date", min_date=numeric_date("2020-01-02"))
    {'strain2'}

    If the requested date column does not exist, we quietly skip this filter.

    >>> sorted(filter_by_min_date(metadata, date_column="missing_column", min_date=numeric_date("2020-01-02")))
    ['strain1', 'strain2']

    """
    strains = set(metadata.index.values)

    # Skip this filter if the date column does not exist.
    if date_column not in metadata.columns:
        return strains

    dates = get_numerical_dates(metadata, date_col=date_column, fmt="%Y-%m-%d")

    filtered = {s for s in strains
                if (dates[s] is not None
                    and (np.isscalar(dates[s]) or all(dates[s]))
                    and np.max(dates[s]) >= min_date)}

    return filtered


def filter_by_max_date(metadata, date_column, max_date) -> FilterFunctionReturn:
    """Filter metadata by maximum date.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    date_column : str
        Column in the dataframe with dates.
    max_date : float
        Maximum date

    Examples
    --------
    >>> from augur.dates import numeric_date
    >>> metadata = pd.DataFrame([{"region": "Africa", "date": "2020-01-01"}, {"region": "Europe", "date": "2020-01-02"}], index=["strain1", "strain2"])
    >>> filter_by_max_date(metadata, date_column="date", max_date=numeric_date("2020-01-01"))
    {'strain1'}

    If the requested date column does not exist, we quietly skip this filter.

    >>> sorted(filter_by_max_date(metadata, date_column="missing_column", max_date=numeric_date("2020-01-01")))
    ['strain1', 'strain2']

    """
    strains = set(metadata.index.values)

    # Skip this filter if the date column does not exist.
    if date_column not in metadata.columns:
        return strains

    dates = get_numerical_dates(metadata, date_col=date_column, fmt="%Y-%m-%d")

    filtered = {s for s in strains
                if (dates[s] is not None
                    and (np.isscalar(dates[s]) or all(dates[s]))
                    and np.min(dates[s]) <= max_date)}

    return filtered


def filter_by_sequence_index(metadata, sequence_index) -> FilterFunctionReturn:
    """Filter metadata by presence of corresponding entries in a given sequence
    index. This filter effectively intersects the strain ids in the metadata and
    sequence index.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    sequence_index : pandas.DataFrame
        Sequence index

    Examples
    --------
    >>> metadata = pd.DataFrame([{"region": "Africa", "date": "2020-01-01"}, {"region": "Europe", "date": "2020-01-02"}], index=["strain1", "strain2"])
    >>> sequence_index = pd.DataFrame([{"strain": "strain1", "ACGT": 28000}]).set_index("strain")
    >>> filter_by_sequence_index(metadata, sequence_index)
    {'strain1'}

    """
    metadata_strains = set(metadata.index.values)
    sequence_index_strains = set(sequence_index.index.values)

    return metadata_strains & sequence_index_strains


def filter_by_sequence_length(metadata, sequence_index, min_length) -> FilterFunctionReturn:
    """Filter metadata by sequence length from a given sequence index.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    sequence_index : pandas.DataFrame
        Sequence index
    min_length : int
        Minimum number of standard nucleotide characters (A, C, G, or T) in each sequence

    Examples
    --------
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


def filter_by_non_nucleotide(metadata, sequence_index) -> FilterFunctionReturn:
    """Filter metadata for strains with invalid nucleotide content.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    sequence_index : pandas.DataFrame
        Sequence index

    Examples
    --------
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


def force_include_strains(metadata, include_file) -> FilterFunctionReturn:
    """Include strains in the given text file from the given metadata.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata indexed by strain name
    include_file : str
        Filename with strain names to include from the given metadata

    Examples
    --------
    >>> import os
    >>> from tempfile import NamedTemporaryFile
    >>> metadata = pd.DataFrame([{"region": "Africa"}, {"region": "Europe"}], index=["strain1", "strain2"])
    >>> with NamedTemporaryFile(delete=False) as include_file:
    ...     characters_written = include_file.write(b'strain1')
    >>> force_include_strains(metadata, include_file.name)
    {'strain1'}
    >>> os.unlink(include_file.name)

    """
    included_strains = read_strains(include_file)
    return set(metadata.index.values) & included_strains


def force_include_where(metadata, include_where) -> FilterFunctionReturn:
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

    Examples
    --------
    >>> metadata = pd.DataFrame([{"region": "Africa"}, {"region": "Europe"}], index=["strain1", "strain2"])
    >>> force_include_where(metadata, "region!=Europe")
    {'strain1'}
    >>> force_include_where(metadata, "region=Europe")
    {'strain2'}
    >>> force_include_where(metadata, "region=europe")
    {'strain2'}

    If the column referenced in the given query does not exist, skip the filter.

    >>> force_include_where(metadata, "missing_column=value")
    set()

    """
    column, op, value = _parse_filter_query(include_where)

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


def construct_filters(args, sequence_index) -> Tuple[List[FilterOption], List[FilterOption]]:
    """Construct lists of filters and inclusion criteria based on user-provided
    arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Command line arguments provided by the user.
    sequence_index : pandas.DataFrame
        Sequence index for the provided arguments.
    """
    exclude_by: List[FilterOption] = []
    include_by: List[FilterOption] = []

    # Force include sequences specified in file(s).
    if args.include:
        # Collect the union of all given strains to include.
        for include_file in args.include:
            include_by.append((
                force_include_strains,
                {
                    "include_file": include_file,
                }
            ))

    # Add sequences with particular metadata attributes.
    if args.include_where:
        for include_where in args.include_where:
            include_by.append((
                force_include_where,
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
                "date_column": METADATA_DATE_COLUMN,
                "ambiguity": args.exclude_ambiguous_dates_by,
            }
        ))

    # Filter by min/max date.
    if args.min_date:
        exclude_by.append((
            filter_by_min_date,
            {
                "min_date": args.min_date,
                "date_column": METADATA_DATE_COLUMN,
            }
        ))
    if args.max_date:
        exclude_by.append((
            filter_by_max_date,
            {
                "max_date": args.max_date,
                "date_column": METADATA_DATE_COLUMN,
            }
        ))

    # Filter by sequence length.
    if args.min_length:
        # Skip VCF files and warn the user that the min length filter does not
        # make sense for VCFs.
        is_vcf = filename_is_vcf(args.sequences)

        if is_vcf: #doesn't make sense for VCF, ignore.
            print_err("WARNING: Cannot use min_length for VCF files. Ignoring...")
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

    if args.group_by:
        # The order in which these are applied later should be broad → specific
        # ambiguity (e.g. year then month), otherwise broad ambiguity will be
        # captured by specific ambiguity.
        if {constants.DATE_YEAR_COLUMN, constants.DATE_MONTH_COLUMN, constants.DATE_WEEK_COLUMN} & set(args.group_by):
            exclude_by.append((
                skip_group_by_with_ambiguous_year,
                {"date_column": METADATA_DATE_COLUMN}
            ))
        if {constants.DATE_MONTH_COLUMN, constants.DATE_WEEK_COLUMN} & set(args.group_by):
            exclude_by.append((
                skip_group_by_with_ambiguous_month,
                {"date_column": METADATA_DATE_COLUMN}
            ))
        if constants.DATE_WEEK_COLUMN in args.group_by:
            exclude_by.append((
                skip_group_by_with_ambiguous_day,
                {"date_column": METADATA_DATE_COLUMN}
            ))

    return exclude_by, include_by


def apply_filters(metadata, exclude_by: List[FilterOption], include_by: List[FilterOption]):
    """Apply a list of filters to exclude or force-include records from the given
    metadata and return the strains to keep, to exclude, and to force include.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata to filter

    Returns
    -------
    set :
        Strains to keep (those that passed all filters)
    list of dict :
        Strains to exclude along with the function that filtered them and the arguments used to run the function.
    list of dict :
        Strains to force-include along with the function that filtered them and the arguments used to run the function.


    For example, filter data by minimum date, but force the include of strains
    from Africa.

    Examples
    --------
    >>> from augur.dates import numeric_date
    >>> metadata = pd.DataFrame([{"region": "Africa", "date": "2020-01-01"}, {"region": "Europe", "date": "2020-10-02"}, {"region": "North America", "date": "2020-01-01"}], index=["strain1", "strain2", "strain3"])
    >>> exclude_by = [(filter_by_min_date, {"date_column": "date", "min_date": numeric_date("2020-04-01")})]
    >>> include_by = [(force_include_where, {"include_where": "region=Africa"})]
    >>> strains_to_keep, strains_to_exclude, strains_to_include = apply_filters(metadata, exclude_by, include_by)
    >>> strains_to_keep
    {'strain2'}
    >>> sorted(strains_to_exclude, key=lambda record: record["strain"])
    [{'strain': 'strain1', 'filter': 'filter_by_min_date', 'kwargs': '[["date_column", "date"], ["min_date", 2020.25]]'}, {'strain': 'strain3', 'filter': 'filter_by_min_date', 'kwargs': '[["date_column", "date"], ["min_date", 2020.25]]'}]
    >>> strains_to_include
    [{'strain': 'strain1', 'filter': 'force_include_where', 'kwargs': '[["include_where", "region=Africa"]]'}]

    We also want to filter by characteristics of the sequence data that we've
    annotated in a sequence index.

    >>> sequence_index = pd.DataFrame([{"strain": "strain1", "A": 7000, "C": 7000, "G": 7000, "T": 7000}, {"strain": "strain2", "A": 6500, "C": 6500, "G": 6500, "T": 6500}, {"strain": "strain3", "A": 1250, "C": 1250, "G": 1250, "T": 1250}]).set_index("strain")
    >>> exclude_by = [(filter_by_sequence_length, {"sequence_index": sequence_index, "min_length": 27000})]
    >>> include_by = [(force_include_where, {"include_where": "region=Europe"})]
    >>> strains_to_keep, strains_to_exclude, strains_to_include = apply_filters(metadata, exclude_by, include_by)
    >>> strains_to_keep
    {'strain1'}
    >>> sorted(strains_to_exclude, key=lambda record: record["strain"])
    [{'strain': 'strain2', 'filter': 'filter_by_sequence_length', 'kwargs': '[["min_length", 27000]]'}, {'strain': 'strain3', 'filter': 'filter_by_sequence_length', 'kwargs': '[["min_length", 27000]]'}]
    >>> strains_to_include
    [{'strain': 'strain2', 'filter': 'force_include_where', 'kwargs': '[["include_where", "region=Europe"]]'}]

    """
    strains_to_keep = set(metadata.index.values)
    strains_to_filter = []
    strains_to_force_include = []
    distinct_strains_to_force_include: Set = set()

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
            include_kwargs_str = _filter_kwargs_to_str(include_kwargs)
            for strain in passed:
                strains_to_force_include.append({
                    "strain": strain,
                    "filter": include_name,
                    "kwargs": include_kwargs_str,
                })

    for filter_function, filter_kwargs in exclude_by:
        # Apply the current function with its given arguments. Each function
        # returns a set of strains that passed the corresponding filter.
        try:
            passed = metadata.pipe(
                filter_function,
                **filter_kwargs,
            )
        except Exception as e:
            if filter_function is filter_by_query:
                if isinstance(e, PandasUndefinedVariableError):
                    raise AugurError(f"Query contains a column that does not exist in metadata.") from e
                raise AugurError(f"Internal Pandas error when applying query:\n\t{e}\nEnsure the syntax is valid per <https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query>.") from e
            else:
                raise

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
            filter_kwargs_str = _filter_kwargs_to_str(filter_kwargs)
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


def _filter_kwargs_to_str(kwargs: FilterFunctionKwargs):
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

    Examples
    --------
    >>> from augur.dates import numeric_date
    >>> from augur.filter.include_exclude_rules import filter_by_sequence_length, filter_by_min_date
    >>> sequence_index = pd.DataFrame([{"strain": "strain1", "ACGT": 28000}, {"strain": "strain2", "ACGT": 26000}, {"strain": "strain3", "ACGT": 5000}]).set_index("strain")
    >>> exclude_by = [(filter_by_sequence_length, {"sequence_index": sequence_index, "min_length": 27000})]
    >>> _filter_kwargs_to_str(exclude_by[0][1])
    '[["min_length", 27000]]'
    >>> exclude_by = [(filter_by_min_date, {"date_column": "date", "min_date": numeric_date("2020-03-01")})]
    >>> _filter_kwargs_to_str(exclude_by[0][1])
    '[["date_column", "date"], ["min_date", 2020.17]]'

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
