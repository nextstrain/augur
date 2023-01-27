import calendar
import datetime
import functools
import re

from .errors import InvalidDate, InvalidYearBounds


def tuple_to_date(year, month, day):
    month = min(month, 12)
    day = min(day, max_day_for_year_month(year, month))

    return datetime.date(year=year, month=month, day=day)


def max_day_for_year_month(year, month):
    return calendar.monthrange(year, month)[1]


def resolve_uncertain_int(uncertain_string, min_or_max):
    """
    Takes a string representation of an integer with uncertain places
    occupied by the character `X`. Returns the minimum or maximum
    possible integer.
    """
    if min_or_max == "min":
        result = int(uncertain_string.replace("X", "0"))
    elif min_or_max == "max":
        result = int(uncertain_string.replace("X", "9"))
    else:
        raise "Tried to resolve an uncertain integer to something other than `min` or `max`."

    if result == 0:
        # A date component cannot be 0. Well, year can, but...
        result = 1

    return result


# This was originally from treetime.utils.ambiguous_date_to_date_range.
class AmbiguousDate:
    """Transforms a date string with uncertainty into the range of possible dates."""

    def __init__(self, uncertain_date, fmt="%Y-%m-%d"):
        self.uncertain_date = uncertain_date
        self.fmt = fmt

        self.assert_only_less_significant_uncertainty()

    def range(self, min_max_year=None):
        """Return the range of possible dates defined by the ambiguous date.

        Impose an upper limit of today's date.
        """
        min_date = tuple_to_date(
            resolve_uncertain_int(self.uncertain_date_components["Y"], "min"),
            resolve_uncertain_int(self.uncertain_date_components["m"], "min"),
            resolve_uncertain_int(self.uncertain_date_components["d"], "min"),
        )

        max_date = tuple_to_date(
            resolve_uncertain_int(self.uncertain_date_components["Y"], "max"),
            resolve_uncertain_int(self.uncertain_date_components["m"], "max"),
            resolve_uncertain_int(self.uncertain_date_components["d"], "max"),
        )

        # Limit dates with ambiguous years to the given bounds.
        if "X" in self.uncertain_date_components["Y"] and min_max_year:
            lower_bound, upper_bound = get_bounds(min_max_year)
            if lower_bound:
                # lower_bound should always be truth-y, but add indentation here for readability.
                if max_date < lower_bound:
                    raise InvalidDate(self.uncertain_date, f"Not possible for date to fall within bounds [{lower_bound}, {upper_bound}]")

                if min_date < lower_bound:
                    min_date = lower_bound

            if upper_bound:
                if upper_bound < min_date:
                    raise InvalidDate(self.uncertain_date, f"Not possible for date to fall within bounds [{lower_bound}, {upper_bound}]")

                if max_date > upper_bound:
                    max_date = upper_bound

        # Limit the min and max dates to be no later than today's date.
        min_date = min(min_date, datetime.date.today())
        max_date = min(max_date, datetime.date.today())

        return (min_date, max_date)

    @property
    @functools.lru_cache()
    def uncertain_date_components(self):
        matches = re.search(self.regex, self.uncertain_date)

        if matches is None:
            raise InvalidDate(self.uncertain_date,
                f"Date does not match format `{self.fmt}`."
            )

        return dict(zip(self.fmt_components, matches.groups()))

    @property
    @functools.lru_cache()
    def fmt_components(self):
        # The `re` module doesn't capture repeated groups, so we'll do it without regexes
        return [component[0] for component in self.fmt.split("%") if len(component) > 0]

    @property
    def regex(self):
        """
        Returns regex defined by the format string.
        Currently only supports %Y, %m, and %d.
        """
        return re.compile(
            "^"
            + self.fmt.replace("%Y", "(....)")
            .replace("%m", "(..?)")
            .replace("%d", "(..?)")
            + "$"
        )

    def assert_only_less_significant_uncertainty(self):
        """
        Raise an exception if a constrained digit appears in a less-significant place
        than an uncertain digit.

        Assuming %Y-%m-%d, these patterns are valid:
            2000-01-01
            2000-01-XX
            2000-XX-XX

        but this is invalid, because month is uncertain but day is constrained:
            2000-XX-01

        These invalid cases are assumed to be unintended use of the tool.
        """
        if "X" in self.uncertain_date_components["Y"]:
            if (
                self.uncertain_date_components["m"] != "XX"
                or self.uncertain_date_components["d"] != "XX"
            ):
                raise InvalidDate(self.uncertain_date,
                    "Year contains uncertainty, so month and day must also be uncertain."
                )
        elif "X" in self.uncertain_date_components["m"]:
            if self.uncertain_date_components["d"] != "XX":
                raise InvalidDate(self.uncertain_date,
                    "Month contains uncertainty, so day must also be uncertain."
                )


def get_bounds(min_max_year):
    """Get exact date bounds based on given years."""
    # This must be an iterable with at least one value.
    assert len(min_max_year) > 0

    if len(min_max_year) > 2:
        raise InvalidYearBounds(f"The year bounds {min_max_year!r} must have only one (lower) or two (lower, upper) bounds.")

    lower_year = int(min_max_year[0])

    if len(min_max_year) == 2:
        upper_year = int(min_max_year[1])
    else:
        upper_year = None

    # Ensure years are properly ordered.
    if lower_year and upper_year and lower_year > upper_year:
        lower_year, upper_year = upper_year, lower_year

    try:
        lower_bound = datetime.date(lower_year, 1, 1)
    except ValueError as error:
        if str(error).startswith("year"):
            raise InvalidYearBounds(f"{lower_year} is not a valid year.") from error
        else:
            raise
    if upper_year:
        try:
            upper_bound = datetime.date(upper_year, 12, 31)
        except ValueError as error:
            if str(error).startswith("year"):
                raise InvalidYearBounds(f"{upper_year} is not a valid year.") from error
            else:
                raise
    else:
        upper_bound = None

    return (lower_bound, upper_bound)
