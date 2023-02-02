import datetime
from typing import List
import treetime.utils
from augur.dates import directives
from augur.dates.directives import directive_is_included
from augur.dates.errors import NotAnExactDate, UnsupportedDateFormat


class ExactDate:
    """Represents a scalar exact date from a date string.

    The date string can be numeric or in a custom format.

    All exact dates can be represented in numeric format, and dates >=
    0001-01-01 can be represented as a datetime.datetime object.
    """
    def __init__(self, date: str, possible_formats: List[str] = None):
        # Original date string.
        self.original_date = date

        # Possible custom formats of the date string in addition to numeric format.
        self.possible_formats = possible_formats
        if not self.possible_formats:
            self.possible_formats = []
        self.validate_formats()

        self.as_numeric = None
        self.parse()

    def validate_formats(self):
        """Validate date formats."""
        for format in self.possible_formats:
            if not directive_is_included(directives.EXACT_DATE, format):
                raise UnsupportedDateFormat()
            # FIXME: check that there is at least one non-directive character in the format string to avoid year-only false positives

    def parse(self):
        """Parse the date string. Raise an error if it does not encode an exact date."""
        # Date strings parsable as floats are exact dates except when also
        # parsable as an integer.
        self.check_year_only()
        if self.parse_numeric():
            return
        if self.parse_custom_format():
            return
        raise NotAnExactDate()

    def check_year_only(self):
        """Raise an error if the date is parsable as an integer.

        This is because an integer date is treated as a year-only ambiguous date."""
        try:
            int(self.original_date)
            raise NotAnExactDate()
        except ValueError:
            pass

    def parse_numeric(self):
        """Try parsing the date directly as a float.
        
        Returns a boolean indicating whether the parsing was successful."""
        try:
            self.as_numeric = float(self.original_date)
            return True
        except ValueError:
            return False

    def parse_custom_format(self):
        """Try parsing the date from one of the possible formats.

        Returns a boolean indicating whether the parsing was successful."""
        for format in self.possible_formats:
            try:
                date = datetime.datetime.strptime(self.original_date, format)
            except ValueError:
                continue
            self.as_numeric = treetime.utils.numeric_date(date)
            return True
        return False

    @property
    def as_datetime(self):
        if self.as_numeric < 1:
            raise ValueError("Date is below the bounds of what can be represented by a datetime object.")
        return treetime.utils.datetime_from_numeric(self.as_numeric)

    def __str__(self):
        try:
            return self.as_datetime.strftime("%Y-%m-%d")
        except ValueError:
            return str(self.as_numeric)

    def __lt__(self, other):
        return self.as_numeric < other.as_numeric

    def __gt__(self, other):
        return self.as_numeric > other.as_numeric
    
    def __eq__(self, other):
        return self.as_numeric == other.as_numeric


def from_datetime(date: datetime.date, format: str):
    return ExactDate(date.strftime(format), possible_formats=[format])
