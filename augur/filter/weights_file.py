from typing import List, Optional

import numpy as np
import pandas as pd

from textwrap import dedent
from augur.errors import AugurError


WEIGHTS_COLUMN = 'weight'
COLUMN_VALUE_FOR_DEFAULT_WEIGHT = 'default'


class InvalidWeightsFile(AugurError):
    def __init__(self, file, error_message):
        super().__init__(f"Bad weights file {file!r}.\n{error_message}")


def read_weights_file(weights_file):
    weights = pd.read_csv(weights_file, delimiter='\t', comment='#', dtype="string")

    if non_numeric_weight_lines := [index + 2 for index in weights[~weights[WEIGHTS_COLUMN].str.lstrip("-").str.isnumeric()].index.tolist()]:
        raise InvalidWeightsFile(weights_file, dedent(f"""\
            Found non-numeric weights on the following lines: {non_numeric_weight_lines}
            {WEIGHTS_COLUMN!r} column must be numeric."""))

    # Cast weights to numeric for calculations
    weights[WEIGHTS_COLUMN] = pd.to_numeric(weights[WEIGHTS_COLUMN])

    if negative_weight_lines := [index + 2 for index in weights[weights[WEIGHTS_COLUMN] < 0].index.tolist()]:
        raise InvalidWeightsFile(weights_file, dedent(f"""\
            Found negative weights on the following lines: {negative_weight_lines}
            {WEIGHTS_COLUMN!r} column must be non-negative."""))

    return weights


def get_weighted_columns(weights_file):
    with open(weights_file) as f:
        has_rows = False
        for row in f:
            has_rows = True
            if row.startswith('#'):
                continue
            columns = row.rstrip().split('\t')
            break
    if not has_rows:
        raise InvalidWeightsFile(weights_file, "File is empty.")
    columns.remove(WEIGHTS_COLUMN)
    return columns


def get_default_weight(weights: pd.DataFrame, weighted_columns: List[str]) -> Optional[np.number]:
    # Match weighted columns with 'default' value. Multiple values can be matched for 2 reasons:
    # 1. Repeated rows following additional permutation with unweighted columns.
    #    This is handled by unique() since the value should be the same.
    # 2. Multiple default rows specified in the weights file.
    #    This is a user error.
    mask = (
        weights[weighted_columns].eq(COLUMN_VALUE_FOR_DEFAULT_WEIGHT).all(axis=1) &
        weights[weighted_columns].notna().all(axis=1)
    )
    default_weight_values = weights.loc[mask, WEIGHTS_COLUMN].unique()  # type: ignore[operator]

    if len(default_weight_values) > 1:
        # TODO: raise InvalidWeightsFile, not AugurError. This function takes
        # the weights DataFrame instead of the filepath, so it does not have the
        # file parameter to InvalidWeightsFile. I didn't want to pass an extra
        # filepath parameter to this function just to have it available for the
        # custom exception class. I also didn't want to pass the filepath
        # parameter. One idea would be to define a custom class to represent a
        # weights file, however this seemed overkill in a quick prototype:
        # <https://github.com/nextstrain/augur/commit/08408580a6beed91e46e93cee28d2ebd552bdbba>
        raise AugurError(f"Multiple default weights were specified: {', '.join(repr(weight) for weight in default_weight_values)}. Only one default weight entry can be accepted.")
    if len(default_weight_values) == 1:
        return default_weight_values[0]
    return None
