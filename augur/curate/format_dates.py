#!/usr/bin/env python3
"""
Standardizes format of date fields of the NDJSON record from stdin to
ISO 8601 date (YYYY-MM-DD) and outputs modified records to stdout.
"""
import argparse
import json
from datetime import datetime
from sys import stderr, stdin, stdout


def format_date(date_string: str, expected_formats: list) -> str:
    """
    Originally from nextstrain/ncov-ingest

    Format *date_string* to ISO 8601 date (YYYY-MM-DD).
    If *date_string* does not match *expected_formats*, return *date_string*.
    If *date_string* is missing the year, return masked date 'XXXX-XX-XX'.
    If *date_string* is an incomplete date (i.e. missing month or day), then
    missing values are masked with 'XX'.

    >>> expected_formats = ['%Y-%m-%d', '%Y-%m-%dT%H:%M:%SZ', '%m-%d']

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
    # Potential directives that datetime accepts that can return the correct year, month, day fields
    # see https://docs.python.org/3.9/library/datetime.html#strftime-and-strptime-format-codes
    #
    # Allows us to check if year/month/day are included in the date format so we
    # know when to mask incomplete dates with 'XX'
    all_field_directives = {'%c', '%x',
        ('%G', '%V', '%A'), ('%G', '%V', '%a'), ('%G', '%V', '%w'), ('%G', '%V', '%u')
    }
    month_and_day_directives = {'%j',
        ('%U', '%A'), ('%U', '%a'), ('%U', '%w'), ('%U', '%u'),
        ('%W', '%A'), ('%W', '%a'), ('%W', '%w'), ('%W', '%u')
    }
    year_directives = {'%y', '%Y'}
    month_directives = {'%b', '%B', '%m'}
    day_directives = {'%d'}

    def directive_is_included(potential_directives: set, date_format: str) -> bool:
        """
        Checks if any of the directives in *potential_directives* is included
        in *date_format* string.

        If an element within *potential_directives* is a tuple, then all directives
        within the tuple must be included in *date_format*.
        """
        return any(
            (
                (isinstance(directive, str) and directive in date_format) or
                (isinstance(directive, tuple) and all(sub_directive in date_format for sub_directive in directive))
            )
            for directive in potential_directives
        )

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

        # If a directive for ALL fields is included in date format,
        # then use all of the parsed field strings
        if (directive_is_included(all_field_directives, date_format)):
            year_string = parsed_year_string
            month_string = parsed_month_string
            day_string = parsed_day_string

        # If not all fields directives are included, then check year
        # directive was included in date format
        elif(directive_is_included(year_directives, date_format)):
            year_string = parsed_year_string

            # Only check for month and day directives if year is included
            # Check if directive for BOTH month and year is included in date format
            if (directive_is_included(month_and_day_directives, date_format)):
                month_string = parsed_month_string
                day_string = parsed_day_string

            # If not directives for BOTH month and day are included, then check
            # month directive was included in date format
            elif(directive_is_included(month_directives, date_format)):
                month_string = parsed_month_string

                # Only check for day directives if month is included
                if(directive_is_included(day_directives, date_format)):
                    day_string = parsed_day_string

        return f"{year_string}-{month_string}-{day_string}"

    if date_string:
        print(
            f"WARNING: Unable to transform date string {date_string!r} because it does not match",
            f"any of the expected formats {expected_formats}.",
            file=stderr
        )

    return date_string


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--date-fields", nargs="+",
        help="List of date field names in the NDJSON record that need to be standardized.")
    parser.add_argument("--expected-date-formats", nargs="+",
        help="Expected date formats that are currently in the provided date fields." +
             "If a date string matches multiple formats, it will be parsed as the first format in the list.")

    args = parser.parse_args()

    expected_formats = args.expected_date_formats

    for record in stdin:
        record = json.loads(record)

        for field in args.date_fields:
            date_string = record.get(field)
            if date_string:
                record[field] = format_date(date_string, expected_formats)

        json.dump(record, stdout, allow_nan=False, indent=None, separators=',:')
        print()
