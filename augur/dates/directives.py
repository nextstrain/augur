from itertools import product
import re


year = {
    # '%y',
    '%Y',
}
month = {
    # '%b',
    # '%B',
    '%m',
}
day = {
    '%d'
}
# month_and_day = {'%j'}
# week = {'%U', '%W'}
# day_of_week = {'%A', '%a', '%w', '%u'}


EXACT_DATE = (
    # Locale's full date representation
    # {('%c',),('%x',)} |

    # Dates with ISO 8601 week dates for year ('%G' is NOT interchangeable with '%Y'), ISO 8601 week ('%V'), and weekdays
    # {('%G', '%V', '%A'),('%G', '%V', '%a'),('%G', '%V', '%w'),('%G', '%V', '%u')} |

    # Dates with year, week, and weekday
    # set(product(year, week, day_of_week)) |

    # Dates with year and day of the year
    # set(product(year, month_and_day)) |

    # Dates with year, month, and day
    set(product(year, month, day))
)


# Set of directives that can be converted to incomplete dates, missing the day
AMBIGUOUS_DAY = set(product(year, month))


# Set of directives that can be converted to incomplete dates, missing the month and day
AMBIGUOUS_MONTH_AND_DAY = set(product(year))



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
