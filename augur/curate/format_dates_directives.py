from itertools import product

year = {'%y', '%Y'}
month = {'%b', '%B', '%m'}
day = {'%d'}
month_and_day = {'%j'}
week = {'%U', '%W'}
day_of_week = {'%A', '%a', '%w', '%u'}

# Set of directives that can be converted to complete date with year, month, and day
YEAR_MONTH_DAY_DIRECTIVES = (
    # Locale's full date representation
    {('%c',),('%x',)} |
    # Dates with ISO 8601 week dates for year ('%G' is NOT interchangeable with '%Y'), ISO 8601 week ('%V'), and weekdays
    {('%G', '%V', '%A'),('%G', '%V', '%a'),('%G', '%V', '%w'),('%G', '%V', '%u')} |
    # Dates with year, week, and weekday
    set(product(year, week, day_of_week)) |
    # Dates with year and day of the year
    set(product(year, month_and_day)) |
    # Dates with year, month, and day
    set(product(year, month, day))
)

# Set of directives that can be converted to incomplete dates, missing the day
YEAR_MONTH_DIRECTIVES = set(product(year, month))

# Set of directives that can be converted to incomplete dates, missing the month and day
YEAR_DIRECTIVES = set(product(year))
