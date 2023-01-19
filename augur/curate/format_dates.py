"""
Format date fields to ISO 8601 dates (YYYY-MM-DD), where incomplete dates
are masked with 'XX' (e.g. 2023 -> 2023-XX-XX).
"""
from datetime import datetime
from augur.io.print import print_err


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("format-dates",
        parents=[parent_subparsers.shared_parser],
        help=__doc__)

    required = parser.add_argument_group(title="REQUIRED")
    required.add_argument("--date-fields", nargs="+",
        help="List of date field names in the record that need to be standardized.")
    required.add_argument("--expected-date-formats", nargs="+",
        help="Expected date formats that are currently in the provided date fields, " +
             "defined by standard format codes as listed at " +
             "https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes. " +
             "If a date string matches multiple formats, it will be parsed as the first format in the list.")
    return parser


def directive_is_included(potential_directives, date_format):
    """
    Checks if any of the directives in *potential_directives* is included
    in *date_format* string.

    If an element within *potential_directives* is a tuple, then all directives
    within the tuple must be included in *date_format*.

    Parameters
    ----------
    potential_directives: set[tuple[str, ...]]
        Set of potential directives to check
    date_format: str
        Date format string to check for directives

    Returns
    -------
    bool:
        Whether the provided *date_format* includes any of the *potential_directives*


    >>> potential_directives = {('%y', '%b', '%d'), ('%y', '%B', '%d'), ('%y', '%m', '%d'),}
    >>> directive_is_included(potential_directives, '%G-%V-%A')
    False
    >>> directive_is_included(potential_directives, '%y-%m')
    False
    >>> directive_is_included(potential_directives, '%y-%m-%d')
    True
    >>> directive_is_included(potential_directives, '%y-%m-%dT%H:%M:%SZ')
    True
    """
    return any(
        all(sub_directive in date_format for sub_directive in directive)
        for directive in potential_directives
    )


def format_date(date_string, expected_formats):
    """
    Format *date_string* to ISO 8601 date (YYYY-MM-DD) by trying to parse it
    as one of the provided *expected_formats*.

    Parameters
    ----------
    date_string: str
        Date string to format
    expected_formats: list[str]
        List of expected formats for the provided date string

    Returns
    -------
    str :
        Formatted date string.
        If *date_string* does not match *expected_formats*, returns original *date_string*.
        If *date_string* is an incomplete date, the date is masked with 'XX'.
        Dates without year will be formatted as 'XXXX-XX-XX', even if month/day are known.
        Dates without month will be formatted as 'YYYY-XX-XX', even if day is known.
        Dates without day will be formatted as 'YYYY-MM-XX'.


    >>> expected_formats = ['%Y', '%Y-%m', '%Y-%m-%d', '%Y-%m-%dT%H:%M:%SZ', '%m-%d']
    >>> format_date("01-01", expected_formats)
    'XXXX-XX-XX'
    >>> format_date("2020", expected_formats)
    '2020-XX-XX'
    >>> format_date("2020-01", expected_formats)
    '2020-01-XX'
    >>> format_date("2020-1-15", expected_formats)
    '2020-01-15'
    >>> format_date("2020-1-1", expected_formats)
    '2020-01-01'
    >>> format_date("2020-01-15", expected_formats)
    '2020-01-15'
    >>> format_date("2020-01-15T00:00:00Z", expected_formats)
    '2020-01-15'
    """
    # Set of directives that can be converted to complete date with year, month, and day
    year_month_day_directives = {
        # Locale's full date representation
        ('%c',),('%x',),
        # Dates with ISO 8601 week dates for year ('%G' is NOT interchangeable with '%Y'), ISO 8601 week ('%V'), and weekdays
        ('%G', '%V', '%A'),('%G', '%V', '%a'),('%G', '%V', '%w'),('%G', '%V', '%u'),
        # Dates with year, week, and weekday
        ('%y', '%U', '%A'), ('%y', '%U', '%a'), ('%y', '%U', '%w'), ('%y', '%U', '%u'),
        ('%y', '%W', '%A'), ('%y', '%W', '%a'), ('%y', '%W', '%w'), ('%y', '%W', '%u'),
        ('%Y', '%U', '%A'), ('%Y', '%U', '%a'), ('%Y', '%U', '%w'), ('%Y', '%U', '%u'),
        ('%Y', '%W', '%A'), ('%Y', '%W', '%a'), ('%Y', '%W', '%w'), ('%Y', '%W', '%u'),
        # Dates with year and day of the year
        ('%y', '%j'), ('%Y', '%j'),
        # Dates with year, month, and day
        ('%y', '%b', '%d'), ('%y', '%B', '%d'), ('%y', '%m', '%d'),
        ('%Y', '%b', '%d'), ('%Y', '%B', '%d'), ('%Y', '%m', '%d'),
    }

    # Set of directives that can be converted to incomplete dates, missing the day
    year_month_directives = {
        ('%y', '%b'), ('%y', '%B'), ('%y', '%m'),
        ('%Y', '%b'), ('%Y', '%B'), ('%Y', '%m'),
    }

    # Set of directives that can be converted to incomplete dates, missing the month and day
    year_directives = {
        ('%y',), ('%Y',)
    }

    for date_format in expected_formats:
        try:
            parsed_date = datetime.strptime(date_string, date_format)
        except ValueError:
            continue

        # Default to date masked as 'XXXX-XX-XX' so we don't return incorrect dates
        year_string = 'XXXX'
        month_string = day_string = 'XX'

        parsed_year_string = str(parsed_date.year)
        parsed_month_string = str(parsed_date.month).zfill(2)
        parsed_day_string = str(parsed_date.day).zfill(2)

        # If directives for all year,month,day fields are included in date_format,
        # then use all of the parsed field strings
        if directive_is_included(year_month_day_directives, date_format):
            year_string = parsed_year_string
            month_string = parsed_month_string
            day_string = parsed_day_string

        # If directives only include year and month are included in date_format,
        # then only use the parsed year and month field strings
        elif directive_is_included(year_month_directives, date_format):
            year_string = parsed_year_string
            month_string = parsed_month_string

        # If directives only include year in date_format, the only use the
        # parsed year field string
        elif directive_is_included(year_directives, date_format):
            year_string = parsed_year_string

        return f"{year_string}-{month_string}-{day_string}"

    if date_string:
        print_err(
            f"WARNING: Unable to transform date string {date_string!r} because it does not match",
            f"any of the expected formats {expected_formats}."
        )

    return date_string


def run(args, records):
    for record in records:
        record = record.copy()

        for field in args.date_fields:
            date_string = record.get(field)
            if date_string:
                record[field] = format_date(date_string, args.expected_date_formats)

        yield record
