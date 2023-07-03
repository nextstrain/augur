from typing import Tuple, Union


# Most logic in Augur does not strictly require a value for numeric dates, so
# they are optional.
ScalarNumericDate = Union[float, None]

# Numeric dates are either exact or defined by a range of possible dates.
NumericDate = Union[
    ScalarNumericDate,
    Tuple[ScalarNumericDate, ScalarNumericDate]
]
