"""
Format date fields to Augur-compatible formats.

If the provided ``--expected-date-formats`` represent incomplete dates then
the incomplete dates are masked with 'XX'. For example, providing
``%Y`` will allow year only dates to be formatted as ``2023-XX-XX``.
"""
import re
from datetime import datetime
from textwrap import dedent

from augur.argparse_ import ExtendOverwriteDefault, SKIP_AUTO_DEFAULT_IN_HELP
from augur.errors import AugurError
from augur.io.print import print_err, indented_list
from augur.types import DataErrorMethod
from augur.utils import first_line
from .format_dates_directives import YEAR_DIRECTIVES, YEAR_MONTH_DIRECTIVES, YEAR_MONTH_DAY_DIRECTIVES


# Builtin date formats that this command should parse
# without additional input from the user.
BUILTIN_DATE_FORMATS = [
    '%Y-%m-%d',
    '%Y-%m-XX',
    '%Y-XX-XX',
    'XXXX-XX-XX',
]
BUILTIN_RANGE_FORMATS = {
    '%Y-%m-%d/%Y-%m-%d': re.compile(
        r'^([0-9]{4}-[0-9]{2}-[0-9]{2})/([0-9]{4}-[0-9]{2}-[0-9]{2})$'
    ),

    # NCBI Datasets can return dates in 'YYYY/YYYY' format.
    # There is no mention of this in Datasets docs, but it is a valid range
    # format for collection dates under submission guidelines.
    # <https://www.ncbi.nlm.nih.gov/WebSub/html/help/collection-date.html>
    # TODO: support all combinations of '<value>/<value>'?
    '%Y/%Y': re.compile(r'^([0-9]{4})/([0-9]{4})$'),

    # NCBI Datasets can return dates in '[YYYY TO YYYY]' format.
    # There is no mention of this in Datasets docs, but it is documented
    # for NCBI Pathogen Detection which may be related.
    # <https://www.ncbi.nlm.nih.gov/pathogens/pathogens_help/#range-searches>
    '[%Y TO %Y]': re.compile(r'^\[([0-9]{4}) TO ([0-9]{4})\]$'),
}


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("format-dates",
        parents=[parent_subparsers.shared_parser],
        help=first_line(__doc__))

    required = parser.add_argument_group(title="REQUIRED")
    required.add_argument("--date-fields", metavar="NAME", nargs="+", action=ExtendOverwriteDefault,
        help="List of date field names in the record that need to be standardized.")

    optional = parser.add_argument_group(title="OPTIONAL")
    optional.add_argument("--expected-date-formats", metavar="FORMAT", nargs="+", action=ExtendOverwriteDefault,
        help=dedent(f"""\
            Custom date formats for values in the provided date fields, defined by standard
            format codes available at
            <https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes>.
            If a value matches multiple formats, it will be parsed using the first match.
            Use 'XX' to match masked parts of the date (e.g. '%%m/XX/%%Y').
            The following formats are builtin and automatically used:
            {", ".join(repr(x).replace("%", "%%") for x in BUILTIN_DATE_FORMATS)}.
            User-provided values are considered after the builtin formats.

            The following date range formats are also builtin and automatically used:
            {", ".join(repr(range_format).replace("%", "%%") for range_format in BUILTIN_RANGE_FORMATS.keys())}.
            Currently there is no option to specify custom date range formats.
            """))
    optional.add_argument("--failure-reporting",
        type=DataErrorMethod.argtype,
        choices=list(DataErrorMethod),
        default=DataErrorMethod.ERROR_FIRST,
        help="How should failed date formatting be reported.")
    optional.add_argument("--no-mask-failure", dest="mask_failure",
        action="store_false",
        help="Do not mask dates with 'XXXX-XX-XX' and return original date string if date formatting failed. " +
             f"(default: False{SKIP_AUTO_DEFAULT_IN_HELP})")

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
    >>> directive_is_included(potential_directives, '%%y-%m-%d')
    False
    >>> directive_is_included(potential_directives, '%y-%m-%d')
    True
    >>> directive_is_included(potential_directives, '%y-%m-%dT%H:%M:%SZ')
    True
    """
    return any(
        all(
            # Exclude escaped directives (e.g. '%%Y' means literal '%Y' not a four digit year)
            bool(re.search(f"(?<!%){re.escape(sub_directive)}", date_format))
            for sub_directive in directive
        )
        for directive in potential_directives
    )


def format_date(date_string, expected_formats):
    """
    Format *date_string* to an Augur-compatible format.

    >>> format_date("2020-01-15", BUILTIN_DATE_FORMATS)
    '2020-01-15'
    >>> format_date("[2001 TO 2002]", BUILTIN_DATE_FORMATS)
    '2001-01-01/2002-12-31'
    """
    if formatted_date := format_to_iso_date(date_string, expected_formats):
        return formatted_date

    if formatted_range := format_to_iso_interval(date_string):
        return formatted_range

    return None


def format_to_iso_date(date_string, expected_formats):
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
    str or None:
        Formatted date string or None if the parsing of the date string failed.
        If *date_string* is an incomplete date, the date is masked with 'XX'.
        Dates without year will be formatted as 'XXXX-XX-XX', even if month/day are known.
        Dates without month will be formatted as 'YYYY-XX-XX', even if day is known.
        Dates without day will be formatted as 'YYYY-MM-XX'.

    >>> format_to_iso_date("", BUILTIN_DATE_FORMATS)
    'XXXX-XX-XX'
    >>> format_to_iso_date("  ", BUILTIN_DATE_FORMATS)
    'XXXX-XX-XX'
    >>> format_to_iso_date("01-01", ['%m-%d'])
    'XXXX-XX-XX'
    >>> format_to_iso_date("2020", ['%Y'])
    '2020-XX-XX'
    >>> format_to_iso_date("2020-01", ['%Y-%m'])
    '2020-01-XX'
    >>> format_to_iso_date("2020-1-15", ['%Y-%m-%d'])
    '2020-01-15'
    >>> format_to_iso_date("2020-1-1", ['%Y-%m-%d'])
    '2020-01-01'
    >>> format_to_iso_date("2020-01-15", BUILTIN_DATE_FORMATS)
    '2020-01-15'
    >>> format_to_iso_date("2020-01-15T00:00:00Z", ['%Y-%m-%dT%H:%M:%SZ'])
    '2020-01-15'
    >>> format_to_iso_date("[2001 TO 2002]", BUILTIN_DATE_FORMATS) is None
    True
    """

    date_string = date_string.strip()
    if date_string == '':
        return 'XXXX-XX-XX'

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
        if directive_is_included(YEAR_MONTH_DAY_DIRECTIVES, date_format):
            year_string = parsed_year_string
            month_string = parsed_month_string
            day_string = parsed_day_string

        # If directives only include year and month are included in date_format,
        # then only use the parsed year and month field strings
        elif directive_is_included(YEAR_MONTH_DIRECTIVES, date_format):
            year_string = parsed_year_string
            month_string = parsed_month_string

        # If directives only include year in date_format, the only use the
        # parsed year field string
        elif directive_is_included(YEAR_DIRECTIVES, date_format):
            year_string = parsed_year_string

        return f"{year_string}-{month_string}-{day_string}"

    return None


def format_to_iso_interval(date_string):
    """
    Format *date_string* to an ISO 8601 interval, compatible with RE_DATE_RANGE
    in augur.dates.

    >>> format_to_iso_interval("[2001 TO 2002]")
    '2001-01-01/2002-12-31'
    >>> format_to_iso_interval("2001-01-01/2002-12-31")
    '2001-01-01/2002-12-31'
    >>> format_to_iso_interval("2001/2002")
    '2001-01-01/2002-12-31'
    """
    if BUILTIN_RANGE_FORMATS['%Y-%m-%d/%Y-%m-%d'].match(date_string):
        return date_string

    if match := BUILTIN_RANGE_FORMATS['%Y/%Y'].match(date_string):
        start_year, end_year = match.groups()
        return f"{start_year}-01-01/{end_year}-12-31"

    if match := BUILTIN_RANGE_FORMATS['[%Y TO %Y]'].match(date_string):
        start_year, end_year = match.groups()
        return f"{start_year}-01-01/{end_year}-12-31"

    # TODO: make BUILTIN_RANGE_FORMATS extendable like --expected-date-formats

    return None


def run(args, records):
    expected_date_formats = BUILTIN_DATE_FORMATS
    if args.expected_date_formats:
        expected_date_formats.extend(args.expected_date_formats)

    failures = []
    failure_reporting = args.failure_reporting
    failure_suggestion = (
        f"Current expected date formats are {expected_date_formats!r}. "
        "This can be updated with --expected-date-formats. "
        f"The following date range formats are also acceptable: {list(BUILTIN_RANGE_FORMATS)!r}. "
        "Currently there is no option to specify custom date range formats."
    )
    for index, record in enumerate(records):
        record = record.copy()
        record_id = index

        for field in args.date_fields:
            date_string = record.get(field)

            if date_string is None:
                raise AugurError(f"Expected date field {field!r} not found in record {record_id!r}.")

            formatted_date_string = format_date(date_string, expected_date_formats)
            if formatted_date_string is None:
                # Mask failed date formatting before processing error methods
                # to ensure failures are masked even when failures are "silent"
                if args.mask_failure:
                    record[field] = "XXXX-XX-XX"

                if failure_reporting is DataErrorMethod.SILENT:
                    continue

                failure_message = f"Unable to format date string {date_string!r} in field {field!r} of record {record_id!r}."
                if failure_reporting is DataErrorMethod.ERROR_FIRST:
                    raise AugurError(dedent(f"""\
                        {failure_message}
                        {failure_suggestion}"""))

                if failure_reporting is DataErrorMethod.WARN:
                    print_err(f"WARNING: {failure_message}")

                # Keep track of failures for final summary
                failures.append((record_id, field, date_string))
            else:
                record[field] = formatted_date_string

        yield record

    if failure_reporting is not DataErrorMethod.SILENT and failures:
        if failure_reporting is DataErrorMethod.ERROR_ALL:
            raise AugurError(dedent(f"""\
                Unable to format dates for the following (record, field, date string):
                {indented_list(map(repr, failures), "                ")}
                {failure_suggestion}"""))

        elif failure_reporting is DataErrorMethod.WARN:
            print_err(f"WARNING: {failure_suggestion}")

        else:
            raise ValueError(f"Encountered unhandled failure reporting method: {failure_reporting!r}")
