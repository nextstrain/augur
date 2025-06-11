"""
Impose lower and/or upper bounds on dates in a column.

Updated values are are formatted as ISO 8601 intervals.
"""
import argparse
import datetime
from textwrap import dedent
from treetime.utils import datestring_from_numeric
from typing import Any, Dict, Iterable, Tuple, Optional, Union

from augur.dates import date_to_numeric, get_numerical_date_from_value
from augur.errors import AugurError
from augur.io.print import print_err, indented_list
from augur.types import DataErrorMethod
from augur.utils import first_line


TODAY = 'today'


RecordType = Dict[str, Any]
     # maybe Dict[str, str]?


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


def run(args: argparse.Namespace, records: Iterable[RecordType]):
    validate_arguments(args)

    failures = []

    for index, input_record in enumerate(records):
        record = input_record.copy()
        try:
            record[args.date_field] = Record(record, index).get_bounded_date(args)
        except DataError as error:
            if args.failure_reporting is DataErrorMethod.SILENT:
                continue
            if args.failure_reporting is DataErrorMethod.ERROR_FIRST:
                raise error
            if args.failure_reporting is DataErrorMethod.WARN:
                print_err(f"WARNING: {error}")
                failures.append(error)
                continue
            if args.failure_reporting is DataErrorMethod.ERROR_ALL:
                failures.append(error)
                continue
            else:
                raise ValueError(f"Encountered unhandled failure reporting method: {args.failure_reporting!r}")
        yield record

    if args.failure_reporting is not DataErrorMethod.SILENT and failures:
        if args.failure_reporting is DataErrorMethod.ERROR_ALL:
            raise AugurError(dedent(f"""\
                Unable to apply bounds. All errors:
                {indented_list(map(str, failures), "                ")}"""))


def validate_arguments(args: argparse.Namespace):
    if not args.lower_bound and not args.upper_bound:
        raise AugurError("At least one of --lower-bound and --upper-bound is required.")


class Record:
    """
    Helper class to wrap a record, its id, and arguments for ease of error handling and small functions.
    """

    def __init__(self, data: RecordType, record_id: int) -> None:
        self.data = data
        self.id = record_id

    def get_bounded_date(self, args: argparse.Namespace) -> str:
        """
        Returns a date string representing the date with bounds applied. bounded interval for the given record.
        """
        start, end = self.convert_date_to_range(args.date_field)

        lower_bound, upper_bound = self.get_bounds(args.lower_bound, args.upper_bound)

        # If any ends are unbounded, return the original date
        if (start == float("-inf") and lower_bound is None) or \
           (end   == float( "inf") and upper_bound is None):
            return self.data[args.date_field]

        # Error if start or end are out of bounds.
        if lower_bound and start < lower_bound and end < lower_bound:
            self.raise_data_error(
                f"{args.date_field!r}={self.data[args.date_field]!r} "
                f"is earlier than the lower bound of "
                f"{args.lower_bound!r}={self.data[args.lower_bound]!r}"
            )
        if upper_bound and start > upper_bound and end > upper_bound:
            self.raise_data_error(
                f"{args.date_field!r}={self.data[args.date_field]!r} "
                f"is later than the upper bound of "
                f"{args.upper_bound!r}={self.data[args.upper_bound]!r}"
            )

        # If the target date overlaps with the bounds, apply the bounds.
        # The start should be no earlier than the lower bound
        # and the end should be no later than the upper bound.
        if lower_bound:
            start = max(start, lower_bound)
        if upper_bound:
            end = min(end, upper_bound)

        # ISO 8601 interval in <start>/<end> format
        return f"{datestring_from_numeric(start)}/{datestring_from_numeric(end)}"

    def convert_date_to_range(self, date_field: str) -> Tuple[float, float]:
        original_date = self.data.get(date_field)

        if original_date is None:
            self.raise_data_error(
                f"Missing date field {date_field!r}."
            )

        date = get_numerical_date_from_value(original_date, fmt="%Y-%m-%d")

        if date == None:
            self.raise_data_error(
                f"Unable to parse value from {date_field!r} as a date: {original_date!r}. "
                "Consider formatting values with augur curate format-dates before applying bounds."
            )

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

    def get_bounds(self, lower_bound_field_or_value: Optional[str], upper_bound_field_or_value: Optional[str]) -> Tuple[Optional[float], Optional[float]]:
        """
        Returns a tuple representing lower and upper bounds.
        """
        lower_bound: Union[float, Tuple[float, float], None] = None
        upper_bound: Union[float, Tuple[float, float], None] = None

        if lower_bound_field_or_value is not None:
            value = self.data.get(lower_bound_field_or_value)
            if value is not None:
                # Input is a field name
                lower_bound = get_numerical_date_from_value(value, fmt="%Y-%m-%d")
                if lower_bound is None:
                    self.raise_data_error(
                        f"Unable to parse value from {lower_bound_field_or_value!r} as a date: {value!r}. "
                        "Consider formatting values with augur curate format-dates before applying bounds."
                    )
            else:
                # Try parsing as a date
                lower_bound = get_numerical_date_from_value(lower_bound_field_or_value, fmt="%Y-%m-%d")
                if lower_bound is None:
                    raise AugurError(f"Expected --lower-bound to be a field name or date, but got {lower_bound_field_or_value!r}.")

        if upper_bound_field_or_value == TODAY:
            if TODAY in self.data:
                raise AugurError(f"{TODAY!r} is ambiguous as it is both an alias to the current date and a field name.")
            upper_bound = date_to_numeric(datetime.date.today())
        elif upper_bound_field_or_value is not None:
            value = self.data.get(upper_bound_field_or_value)
            if value is not None:
                # Input is a field name
                upper_bound = get_numerical_date_from_value(value, fmt="%Y-%m-%d")
                if upper_bound is None:
                    self.raise_data_error(
                        f"Unable to parse value from {upper_bound_field_or_value!r} as a date: {value!r}. "
                        "Consider formatting values with augur curate format-dates before applying bounds."
                    )
            else:
                # Try parsing as a date
                upper_bound = get_numerical_date_from_value(upper_bound_field_or_value, fmt="%Y-%m-%d")
                if upper_bound is None:
                    raise AugurError(f"Expected --upper-bound to be a field name or date, but got {upper_bound_field_or_value!r}.")

        # Resolve ranges to single values
        if isinstance(lower_bound, tuple):
            lower_bound = lower_bound[0]
        if isinstance(upper_bound, tuple):
            upper_bound = upper_bound[1]

        return lower_bound, upper_bound

    def raise_data_error(self, message: str) -> None:
        raise DataError(self.id, message)


class DataError(AugurError):
    def __init__(self, record_id: int, message: str):
        self.record_id = record_id
        self.message = message

    def __str__(self):
        return f"[record {self.record_id}] {self.message}"
