"""
Impose lower and/or upper bounds on dates in a column.

Updated values are are formatted as ISO 8601 intervals.
"""
from textwrap import dedent
from treetime.utils import datestring_from_numeric

from augur.dates import get_numerical_date_from_value
from augur.errors import AugurError
from augur.io.print import print_err, indented_list
from augur.types import DataErrorMethod
from augur.utils import first_line


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
    # FIXME: support constants
    optional.add_argument("--lower-bound", metavar="NAME",
        help="Name of an existing date field to use as the lower bound for --date-field (i.e. minimum)")
    optional.add_argument("--upper-bound", metavar="NAME",
        help="Name of an existing date field to use as the upper bound for --date-field (i.e. maximum)")
    optional.add_argument("--failure-reporting",
        type=DataErrorMethod.argtype,
        choices=list(DataErrorMethod),
        default=DataErrorMethod.ERROR_FIRST,
        help="How should failed date formatting be reported.")

    return parser


def validate_arguments(args):
    if not args.lower_bound and not args.upper_bound:
        raise AugurError("At least one of --lower-bound and --upper-bound is required.")


def run(args, records):
    validate_arguments(args)

    failure_reporting = args.failure_reporting

    failures = []
    failure_suggestion = ...
    if args.lower_bound and args.upper_bound:
        failure_suggestion = (
            f"Current bounds for incomplete dates in {args.date_field!r} "
            f"are defined by a minimum from {args.lower_bound!r} "
            f"and a maximum from {args.upper_bound!r}. "
            "These can be updated with --lower-bound and --upper-bound."
        )
    elif args.lower_bound:
        failure_suggestion = (
            f"Current bounds for incomplete dates in {args.date_field!r} "
            f"are defined by a minimum from {args.lower_bound!r}. "
            "This can be updated with --lower-bound."
        )
    elif args.upper_bound:
        failure_suggestion = (
            f"Current bounds for incomplete dates in {args.date_field!r} "
            f"are defined by a maximum from {args.upper_bound!r}. "
            "This can be updated with --upper-bound."
        )


    def apply_bounds(record, record_id):
        record = record.copy()

        # Initialize start/end from original date
        original_date = get_numerical_date_from_value(record[args.date_field], fmt="%Y-%m-%d")
        start, end = ..., ...
        if isinstance(original_date, tuple):
            start, end = original_date
        elif isinstance(original_date, float) or isinstance(original_date, int):
            start = original_date
            end = original_date
        else:
            start = float('-inf')
            end = float('inf')

        lower_bound, upper_bound = None, None

        # Get bounds
        if args.lower_bound:
            lower_bound = get_numerical_date_from_value(record[args.lower_bound], fmt="%Y-%m-%d")
            if isinstance(lower_bound, tuple):
                lower_bound = lower_bound[0]
        if args.upper_bound:
            upper_bound = get_numerical_date_from_value(record[args.upper_bound], fmt="%Y-%m-%d")
            if isinstance(upper_bound, tuple):
                upper_bound = upper_bound[1]

        # If any infinities are unbounded, don't modify the date
        if (start == float("-inf") and lower_bound is None) or \
           (end   == float( "inf") and upper_bound is None):
            return record

        # Error if the target date does not overlap with the bounds
        target_out_of_bounds = False

        # Check lower bound
        if lower_bound and start < lower_bound and end < lower_bound:
            failure_message = (
                f"{args.date_field!r}={record[args.date_field]!r} "
                f"is earlier than the lower bound of "
                f"{args.lower_bound!r}={record[args.lower_bound]!r}"
            )
            if failure_reporting is DataErrorMethod.SILENT:
                pass
            elif failure_reporting is DataErrorMethod.ERROR_FIRST:
                raise AugurError(failure_message)
            elif failure_reporting is DataErrorMethod.WARN:
                print_err(f"WARNING: {failure_message}")
                failures.append((record_id, args.date_field, record[args.date_field], lower_bound, upper_bound))
            elif failure_reporting is DataErrorMethod.ERROR_ALL:
                failures.append((record_id, args.date_field, record[args.date_field], lower_bound, upper_bound))
            target_out_of_bounds = True

        # Check upper bound
        if upper_bound and start > upper_bound and end > upper_bound:
            failure_message = (
                f"{args.date_field!r}={record[args.date_field]!r} "
                f"is later than the upper bound of "
                f"{args.upper_bound!r}={record[args.upper_bound]!r}"
            )
            if failure_reporting is DataErrorMethod.SILENT:
                pass
            elif failure_reporting is DataErrorMethod.ERROR_FIRST:
                raise AugurError(failure_message)
            elif failure_reporting is DataErrorMethod.WARN:
                print_err(f"WARNING: {failure_message}")
                failures.append((record_id, args.date_field, record[args.date_field], lower_bound, upper_bound))
            elif failure_reporting is DataErrorMethod.ERROR_ALL:
                failures.append((record_id, args.date_field, record[args.date_field], lower_bound, upper_bound))
            target_out_of_bounds = True

        # FIXME: factor this out into a separate unit-test-able function
        # If the target date overlaps with the bounds, apply the bounds
        if not target_out_of_bounds:
            # The start should be no earlier than the lower bound
            # and the end should be no later than the upper bound
            if lower_bound:
                start = max(start, lower_bound)
            if upper_bound:
                end = min(end, upper_bound)

            # ISO 8601 interval in <start>/<end> format
            record[args.date_field] = (
                f"{datestring_from_numeric(start)}/{datestring_from_numeric(end)}"
            )

        return record

    for index, record in enumerate(records):
        date_string = record.get(args.date_field)
        if date_string is None:
            raise AugurError(f"Expected date field {args.date_field!r} not found in record {record_id!r}.")

        record = apply_bounds(record, index)

        yield record

    if failure_reporting is not DataErrorMethod.SILENT and failures:
        if failure_reporting is DataErrorMethod.ERROR_ALL:
            raise AugurError(dedent(f"""\
                Unable to apply bounds for the following (record, field, formatted date string, lower bound, upper bound):
                {indented_list(map(repr, failures), "                ")}
                {failure_suggestion}"""))

        elif failure_reporting is DataErrorMethod.WARN:
            print_err(f"WARNING: {failure_suggestion}")

        else:
            raise ValueError(f"Encountered unhandled failure reporting method: {failure_reporting!r}")
