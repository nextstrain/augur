import datetime
import re
from typing import List

from .ambiguous_date import AmbiguousDate
from .exact_date import ExactDate
from .errors import InvalidDate, NotAScalarAmbiguousDate, NotAnExactDate


RE_DATE_RANGE = re.compile(r'^\[(.*),(.*)\]$')
DEFAULT_FORMATS = ['%Y-%m-%d', '%Y-%m', '%Y']


class AugurDate:
    """Represents a date for usage in Augur.

    This date can be:
    1. A scalar exact date
    2. A pair of min and max scalar exact dates
    3. A scalar incomplete/ambiguous date

    All Augur dates can be represented in a range format of "[min,max]", with
    min and max being inclusive exact date bounds.
    """
    def __init__(self, date: str, possible_formats: List[str] = None):
        # Original date string.
        self.original_date = date

        # Possible scalar date formats.
        self.possible_formats = possible_formats
        if not self.possible_formats:
            self.possible_formats = DEFAULT_FORMATS

        self.format = None
        self.min = None
        self.max = None

        self.parse()

    @property
    def is_null(self):
        """Represents whether the date is parsed as a null/empty value.

        This happens when a date string is parsed without errors but no date
        range can be determined.
        """
        # I chose to create a separate property instead of implementing
        # self.__bool__() to make it easier to track where this is used.
        return self.min is None and self.max is None

    def parse(self):
        if self.parse_exact():
            return
        if self.parse_range():
            return
        if self.parse_negative_year():
            return

        # At this point, the date can only be valid if it's an ambiguous date in
        # one of the possible formats.
        try:
            if self.parse_ambiguous_scalar():
                return
        except NotAScalarAmbiguousDate:
            return

    def range(self):
        return (self.min, self.max)

    def parse_exact(self):
        """Try parsing the date string as an exact date.

        Returns a boolean indicating whether the parsing was successful."""
        try:
            exact_date = ExactDate(self.original_date, self.possible_formats)
            self.min = exact_date
            self.max = exact_date
            return True
        except NotAnExactDate:
            return False

    def parse_range(self):
        """Try parsing the date string as the range format that is returned by this class.

        Returns a boolean indicating whether the parsing was successful."""
        match = RE_DATE_RANGE.match(self.original_date)
        if match:
            matched_min, matched_max = match.groups()
            try:
                self.min = ExactDate(matched_min, self.possible_formats)
                self.max = ExactDate(matched_max, self.possible_formats)
            except NotAnExactDate:
                raise InvalidDate(self.original_date, "Range provided but at least one end is not a valid exact date.")

            if self.min > self.max:
                raise InvalidDate(self.original_date, "Range must be in chronological order.")

            return True

        return False

    def parse_negative_year(self):
        """Try parsing the date string as a negative year-only ambiguous date.

        Positive year-only scalar dates are supported as ambiguous dates
        representing the range of all days in that particular year. Similar
        support should be added for negative years, though it must be done
        separately since negative dates cannot be parsed by a custom datetime
        format.

        Returns a boolean indicating whether the parsing was successful."""

        assert "." not in self.original_date

        try:
            year = int(self.original_date)
        except ValueError:
            return False

        if year < 0:
            self.min = ExactDate(f"{year}.9999")
            self.max = ExactDate(f"{year}.0")
            return True

        return False

    def parse_ambiguous_scalar(self):
        """Try parsing the date string as a scalar ambiguous date.

        Returns a boolean indicating whether the parsing was successful."""
        incomplete_format = self.find_format_of_ambiguous_date(self.original_date)
        # FIXME: somehow get full format from incomplete format?
        ambiguous_date = AmbiguousDate(self.original_date, incomplete_format)
        # FIXME: handle min_max_year
        self.min, self.max = ambiguous_date.range()
        return True

    def find_format_of_ambiguous_date(self, date_str):
        # Remove any ambiguity to determine format.
        # A replacement value of 1 results in the least amount of possible errors to be handled.
        # Note: A month with ambiguous
        date = date_str.replace('X', '1')

        for date_format in self.possible_formats:
            try:
                datetime.datetime.strptime(date, date_format)
                return date_format
            except ValueError as error:
                if 'day is out of range for month' in str(error):
                    # This can be an artifact of ambiguity or user value.
                    # These can be ignored for the purposes of determining format.
                    return date_format

        # Raise an error if the date matches none of the possible formats.
        raise NotAScalarAmbiguousDate()

    def __str__(self):
        """Return date in [min,max] format."""
        return f"[{self.min},{self.max}]"
