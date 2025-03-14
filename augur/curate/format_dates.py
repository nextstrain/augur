"""
Format date fields to ISO 8601-like formats.

This command provides two main functionalities:

1. Normalize date formats

   When date fields are specified with ``--date-fields``, values are normalized
   to complete or partial ISO 8601 date strings. Values matching only year or
   year-month formats are masked with 'XX' for missing components (e.g., ``%Y``
   → ``2023-XX-XX``). Use ``--expected-date-formats`` to provide custom date
   formats for values in the provided date fields.

2. Apply interval bounds

   When a target date field is specified with ``--target-date-field``,
   incomplete dates in that field will be formatted as ISO 8601 intervals using
   lower and/or upper bounds from
   ``--target-date-field-min``/``--target-date-field-max``.
"""
import re
from datetime import datetime
from textwrap import dedent
from treetime.utils import datestring_from_numeric

from augur.argparse_ import ExtendOverwriteDefault, SKIP_AUTO_DEFAULT_IN_HELP
from augur.dates import get_numerical_date_from_value
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


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("format-dates",
        parents=[parent_subparsers.shared_parser],
        help=first_line(__doc__))

    optional = parser.add_argument_group(title="OPTIONAL")
    optional.add_argument("--date-fields", metavar="NAME", nargs="+", action=ExtendOverwriteDefault,
        help="List of date field names in the record that need to be standardized.")
    optional.add_argument("--expected-date-formats", metavar="FORMAT", nargs="+", action=ExtendOverwriteDefault,
        help=dedent(f"""\
            Custom date formats for values in the provided date fields, defined by standard
            format codes available at
            <https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes>.
            If a value matches multiple formats, it will be parsed using the first match.
            Use 'XX' to match masked parts of the date (e.g. '%%m/XX/%%Y').
            The following formats are builtin and automatically used:
            {", ".join(repr(x).replace("%", "%%") for x in BUILTIN_DATE_FORMATS)}.
            User-provided values are considered after the builtin formats."""))
    optional.add_argument("--failure-reporting",
        type=DataErrorMethod.argtype,
        choices=list(DataErrorMethod),
        default=DataErrorMethod.ERROR_FIRST,
        help="How should failed date formatting be reported.")
    optional.add_argument("--no-mask-failure", dest="mask_failure",
        action="store_false",
        help="Do not mask dates with 'XXXX-XX-XX' and return original date string if date formatting failed. " +
             f"(default: False{SKIP_AUTO_DEFAULT_IN_HELP})")
    optional.add_argument("--target-date-field", metavar="NAME",
        help=dedent("""\
            Name of an existing date field to apply bounds to. Incomplete values
            for this field will be formatted as an interval using bounds
            provided by --target-date-field-min and/or --target-date-field-max."""))
    optional.add_argument("--target-date-field-min", metavar="NAME",
        help="Name of an existing date field to use as the lower bound for --target-date-field (i.e. minimum)")
    optional.add_argument("--target-date-field-max", metavar="NAME",
        help="Name of an existing date field to use as the upper bound for --target-date-field (i.e. maximum)")

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


    >>> expected_formats = ['%Y', '%Y-%m', '%Y-%m-%d', '%Y-%m-%dT%H:%M:%SZ', '%m-%d']
    >>> format_date("", expected_formats)
    'XXXX-XX-XX'
    >>> format_date("  ", expected_formats)
    'XXXX-XX-XX'
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


def validate_arguments(args):
    if args.target_date_field and not args.target_date_field_min and not args.target_date_field_max:
        raise AugurError("--target-date-field requires at least one of --target-date-field-min, --target-date-field-max.")
    if args.target_date_field_min and not args.target_date_field:
        raise AugurError("--target-date-field-min requires --target-date-field.")
    if args.target_date_field_max and not args.target_date_field:
        raise AugurError("--target-date-field-max requires --target-date-field.")


def run(args, records):
    validate_arguments(args)

    expected_date_formats = BUILTIN_DATE_FORMATS
    if args.expected_date_formats:
        expected_date_formats.extend(args.expected_date_formats)

    failure_reporting = args.failure_reporting

    format_failures = []
    format_failure_suggestion = (
        f"Current expected date formats are {expected_date_formats!r}. " +
        "This can be updated with --expected-date-formats."
    )

    bounds_failures = []
    bounds_failure_suggestion = ...
    if args.target_date_field_min and args.target_date_field_max:
        bounds_failure_suggestion = (
            f"Current bounds for incomplete dates in {args.target_date_field!r} "
            f"are defined by a minimum from {args.target_date_field_min!r} "
            f"and a maximum from {args.target_date_field_max!r}. "
            "These can be updated with --target-date-field-min and --target-date-field-max."
        )
    elif args.target_date_field_min:
        bounds_failure_suggestion = (
            f"Current bounds for incomplete dates in {args.target_date_field!r} "
            f"are defined by a minimum from {args.target_date_field_min!r}. "
            "This can be updated with --target-date-field-min."
        )
    elif args.target_date_field_max:
        bounds_failure_suggestion = (
            f"Current bounds for incomplete dates in {args.target_date_field!r} "
            f"are defined by a maximum from {args.target_date_field_max!r}. "
            "This can be updated with --target-date-field-max."
        )

    def normalize_format(record, record_id):
        record = record.copy()

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
                        {format_failure_suggestion}"""))

                if failure_reporting is DataErrorMethod.WARN:
                    print_err(f"WARNING: {failure_message}")

                # Keep track of failures for final summary
                format_failures.append((record_id, field, date_string))
            else:
                record[field] = formatted_date_string

        return record

    def apply_bounds(record, record_id):
        record = record.copy()

        original_date = get_numerical_date_from_value(record[args.target_date_field], fmt="%Y-%m-%d")

        # Keep exact dates as-is
        # Maybe worth adding a warning if the date is out of bounds?
        if not isinstance(original_date, tuple):
            return record

        start, end = original_date
        lower_bound, upper_bound = None, None

        # Get bounds
        if args.target_date_field_min:
            lower_bound = get_numerical_date_from_value(record[args.target_date_field_min], fmt="%Y-%m-%d")
            if isinstance(lower_bound, tuple):
                lower_bound = lower_bound[0]
        if args.target_date_field_max:
            upper_bound = get_numerical_date_from_value(record[args.target_date_field_max], fmt="%Y-%m-%d")
            if isinstance(upper_bound, tuple):
                upper_bound = upper_bound[1]

        # Error if the target date does not overlap with the bounds
        target_out_of_bounds = False

        # Check lower bound
        if lower_bound and start < lower_bound and end < lower_bound:
            failure_message = (
                f"{args.target_date_field!r}={record[args.target_date_field]!r} "
                f"is earlier than the lower bound of "
                f"{args.target_date_field_min!r}={record[args.target_date_field_min]!r}"
            )
            if failure_reporting is DataErrorMethod.SILENT:
                pass
            elif failure_reporting is DataErrorMethod.ERROR_FIRST:
                raise AugurError(failure_message)
            elif failure_reporting is DataErrorMethod.WARN:
                print_err(f"WARNING: {failure_message}")
                bounds_failures.append((record_id, args.target_date_field, record[args.target_date_field], lower_bound, upper_bound))
            elif failure_reporting is DataErrorMethod.ERROR_ALL:
                bounds_failures.append((record_id, args.target_date_field, record[args.target_date_field], lower_bound, upper_bound))
            target_out_of_bounds = True

        # Check upper bound
        if upper_bound and start > upper_bound and end > upper_bound:
            failure_message = (
                f"{args.target_date_field!r}={record[args.target_date_field]!r} "
                f"is later than the upper bound of "
                f"{args.target_date_field_max!r}={record[args.target_date_field_max]!r}"
            )
            if failure_reporting is DataErrorMethod.SILENT:
                pass
            elif failure_reporting is DataErrorMethod.ERROR_FIRST:
                raise AugurError(failure_message)
            elif failure_reporting is DataErrorMethod.WARN:
                print_err(f"WARNING: {failure_message}")
                bounds_failures.append((record_id, args.target_date_field, record[args.target_date_field], lower_bound, upper_bound))
            elif failure_reporting is DataErrorMethod.ERROR_ALL:
                bounds_failures.append((record_id, args.target_date_field, record[args.target_date_field], lower_bound, upper_bound))
            target_out_of_bounds = True

        # If the target date overlaps with the bounds, apply the bounds
        if not target_out_of_bounds:
            # The start should be no earlier than the lower bound
            # and the end should be no later than the upper bound
            if lower_bound:
                start = max(start, lower_bound)
            if upper_bound:
                end = min(end, upper_bound)

            # ISO 8601 interval in <start>/<end> format
            record[args.target_date_field] = (
                f"{datestring_from_numeric(start)}/{datestring_from_numeric(end)}"
            )

        return record

    for index, record in enumerate(records):
        if args.date_fields:
            record = normalize_format(record, index)

            if failure_reporting is DataErrorMethod.WARN and format_failures:
                print_err(f"WARNING: {format_failure_suggestion}")

        # Apply bounds after normalizing format so that any existing ambiguity is in a resolvable format
        if args.target_date_field:
            record = apply_bounds(record, index)

            if failure_reporting is DataErrorMethod.WARN and bounds_failures:
                print_err(f"WARNING: {bounds_failure_suggestion}")

        yield record

    if failure_reporting is DataErrorMethod.ERROR_ALL and (format_failures or bounds_failures):
        failure_message = ""
        if format_failures:
            failure_message += dedent(f"""\
                Unable to normalize format for the following (record, field, date string):
                {indented_list(map(repr, format_failures), "                ")}
                {format_failure_suggestion}""")
        if bounds_failures:
            if failure_message:
                failure_message += "\n"
            failure_message += dedent(f"""\
                Unable to apply bounds for the following (record, field, formatted date string, lower bound, upper bound):
                {indented_list(map(repr, bounds_failures), "                ")}
                {bounds_failure_suggestion}""")
        raise AugurError(failure_message)
