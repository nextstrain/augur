import calendar
import datetime
import functools
import re


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


class DateDisambiguator:
    """Transforms a date string with uncertainty into the range of possible dates."""

    def __init__(self, uncertain_date, fmt="%Y-%m-%d", min_max_year=None):
        self.uncertain_date = uncertain_date
        self.fmt = fmt
        self.min_max_year = min_max_year

        self.assert_only_less_significant_uncertainty()

    def range(self):
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
        max_date = min(max_date, datetime.date.today())

        return (min_date, max_date)

    @property
    @functools.lru_cache()
    def uncertain_date_components(self):
        matches = re.search(self.regex, self.uncertain_date)

        if matches is None:
            raise ValueError(
                f"Malformed uncertain date `{self.uncertain_date}` for format `{self.fmt}`"
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
                raise ValueError(
                    "Invalid date: Year contains uncertainty, so month and day must also be uncertain."
                )
        elif "X" in self.uncertain_date_components["m"]:
            if self.uncertain_date_components["d"] != "XX":
                raise ValueError(
                    "Invalid date: Month contains uncertainty, so day must also be uncertain."
                )
