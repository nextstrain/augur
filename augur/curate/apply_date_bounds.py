"""
Impose lower and/or upper bounds on dates in a column.

Updated values are are formatted as ISO 8601 intervals.
"""
from textwrap import dedent
from treetime.utils import datestring_from_numeric
from typing import Tuple, Optional

from augur.dates import date_to_numeric, get_numerical_date_from_value
from augur.errors import AugurError
from augur.io.print import print_err, indented_list
from augur.types import DataErrorMethod
from augur.utils import first_line


TODAY = 'today'


class DataError(AugurError):
    def __init__(self, message: str, record_id: int, original_date: str,
                 lower_bound: Optional[float], upper_bound: Optional[float]):
        super().__init__(message)
        self.record_id = record_id
        self.original_date = original_date
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("apply-date-bounds",
        parents=[parent_subparsers.shared_parser],
        help=first_line(__doc__))

    required = parser.add_argument_group(title="REQUIRED")
    required.add_argument("--date-field", metavar="NAME",
        help=dedent("""\
            Name of an existing date field to apply bounds to. Values will be
            formatted as an interval using bounds provided by --lower-bound
            and/or --upper-bound."""))

    optional = parser.add_argument_group(title="OPTIONAL")
    optional.add_argument("--lower-bound", metavar="NAME | DATE",
        help=dedent("""\
            Name of an existing date field or date to use as the lower bound for
            --date-field (i.e. minimum)."""))
    optional.add_argument("--upper-bound", metavar=f"NAME | DATE | {TODAY!r}",
        help=dedent(f"""\
            Name of an existing date field or date to use as the upper bound for
            --date-field (i.e. maximum). Use {TODAY!r} to set the current date as
            the upper bound."""))
    optional.add_argument("--failure-reporting",
        type=DataErrorMethod.argtype,
        choices=list(DataErrorMethod),
        default=DataErrorMethod.ERROR_FIRST,
        help="How should failed date formatting be reported.")

    return parser


def validate_arguments(args):
    if not args.lower_bound and not args.upper_bound:
        raise AugurError("At least one of --lower-bound and --upper-bound is required.")


def get_bounds(args, record) -> Tuple[Optional[float], Optional[float]]:
    """
    Returns a tuple representing lower and upper bounds.
    """
    lower_bound, upper_bound = None, None

    if args.lower_bound:
        if value := record.get(args.lower_bound):
            # Valid field name
            lower_bound = get_numerical_date_from_value(value, fmt="%Y-%m-%d")
            if lower_bound is None:
                raise AugurError(f"Unable to parse value from {args.lower_bound!r} as a date.")
        else:
            # Try parsing as date
            lower_bound = get_numerical_date_from_value(args.lower_bound, fmt="%Y-%m-%d")
            if lower_bound is None:
                raise AugurError(f"Expected --lower-bound to be a field name or date, but got {args.lower_bound!r}.")

    if args.upper_bound:
        if args.upper_bound == TODAY:
            if TODAY in record:
                raise AugurError(f"{TODAY!r} is ambiguous as it is both an alias to the current date and a field name.")
            upper_bound = date_to_numeric(datetime.datetime.today())
        elif value := record.get(args.upper_bound):
            # Valid field name
            upper_bound = get_numerical_date_from_value(value, fmt="%Y-%m-%d")
            if upper_bound is None:
                raise AugurError(f"Unable to parse value from {args.upper_bound!r} as a date.")
        else:
            # Try parsing as date
            upper_bound = get_numerical_date_from_value(args.upper_bound, fmt="%Y-%m-%d")
            if upper_bound is None:
                raise AugurError(f"Expected --upper-bound to be a field name or date, but got {args.upper_bound!r}.")

    # Resolve ranges to single values
    if isinstance(lower_bound, tuple):
        lower_bound = lower_bound[0]
    if isinstance(upper_bound, tuple):
        upper_bound = upper_bound[1]

    return lower_bound, upper_bound


def convert_to_range(date) -> Tuple[Optional[float], Optional[float]]:
    """
    Convert a date to a range (start, end).
    """
    date = get_numerical_date_from_value(date, fmt="%Y-%m-%d")

    start, end = float('-inf'), float('inf')

    if isinstance(date, tuple):
        start, end = date
    elif isinstance(date, float):
        start = date
        end = date
    elif isinstance(date, int):
        start = float(date)
        end = float(date)

    return start, end


def run(args, records):
    validate_arguments(args)

    failure_reporting = args.failure_reporting

    failures = []
    failure_suggestion = ...
    if args.lower_bound and args.upper_bound:
        failure_suggestion = (
            f"The lower bound is {args.lower_bound!r} "
            f"and the upper bound is {args.upper_bound!r}. "
            "These can be updated with --lower-bound and --upper-bound."
        )
    elif args.lower_bound:
        failure_suggestion = (
            f"The lower bound is {args.lower_bound!r}. "
            "This can be updated with --lower-bound."
        )
    elif args.upper_bound:
        failure_suggestion = (
            f"The upper bound is {args.upper_bound!r}. "
            "This can be updated with --upper-bound."
        )

    def get_bounded_date(record, record_id):
        """
        Returns a bounded interval for the given record.
        """
        record = record.copy()

        original_date = record.get(args.date_field)
        if original_date is None:
            raise AugurError(f"Expected date field {args.date_field!r} not found in record {record_id!r}.")

        start, end = convert_to_range(original_date)
        # FIXME: suggest augur curate format-dates when dates are not parseable?

        lower_bound, upper_bound = get_bounds(args, record)

        # If any ends are unbounded, leave the date unchanged
        if (start == float("-inf") and lower_bound is None) or \
           (end   == float( "inf") and upper_bound is None):
            return original_date

        # Check lower bound
        if lower_bound and start < lower_bound and end < lower_bound:
            failure_message = (
                f"{args.date_field!r}={record[args.date_field]!r} "
                f"is earlier than the lower bound of "
                f"{args.lower_bound!r}={record[args.lower_bound]!r}"
            )
            raise DataError(
                failure_message,
                record_id,
                record[args.date_field],
                lower_bound,
                upper_bound
            )

        # Check upper bound
        if upper_bound and start > upper_bound and end > upper_bound:
            failure_message = (
                f"{args.date_field!r}={record[args.date_field]!r} "
                f"is later than the upper bound of "
                f"{args.upper_bound!r}={record[args.upper_bound]!r}"
            )
            raise DataError(
                failure_message,
                record_id,
                record[args.date_field],
                lower_bound,
                upper_bound
            )

        # If the target date overlaps with the bounds, apply the bounds
        # The start should be no earlier than the lower bound
        # and the end should be no later than the upper bound
        if lower_bound:
            start = max(start, lower_bound)
        if upper_bound:
            end = min(end, upper_bound)

        # ISO 8601 interval in <start>/<end> format
        return f"{datestring_from_numeric(start)}/{datestring_from_numeric(end)}"

    for index, record in enumerate(records):
        try:
            record[args.date_field] = get_bounded_date(record, index)
        except DataError as error:
            if failure_reporting is DataErrorMethod.SILENT:
                continue
            if failure_reporting is DataErrorMethod.ERROR_FIRST:
                raise error
            if failure_reporting is DataErrorMethod.WARN:
                print_err(f"WARNING: {error}")
                failures.append((error.record_id, error.original_date, error.lower_bound, error.upper_bound))
                continue
            if failure_reporting is DataErrorMethod.ERROR_ALL:
                failures.append((error.record_id, error.original_date, error.lower_bound, error.upper_bound))
                continue
            else:
                raise ValueError(f"Encountered unhandled failure reporting method: {failure_reporting!r}")
        yield record

    if failure_reporting is not DataErrorMethod.SILENT and failures:
        if failure_reporting is DataErrorMethod.ERROR_ALL:
            raise AugurError(dedent(f"""\
                Unable to apply bounds for the following (record, date, lower bound, upper bound):
                {indented_list(map(repr, failures), "                ")}
                {failure_suggestion}"""))

        elif failure_reporting is DataErrorMethod.WARN:
            print_err(f"WARNING: {failure_suggestion}")
