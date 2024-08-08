import pandas as pd
from textwrap import dedent
from augur.errors import AugurError


WEIGHTS_COLUMN = 'weight'


class InvalidWeightsFile(AugurError):
    def __init__(self, file, error_message):
        super().__init__(f"Bad weights file {file!r}.\n{error_message}")


def read_weights_file(weights_file):
    weights = pd.read_csv(weights_file, delimiter='\t', comment='#')

    if not pd.api.types.is_numeric_dtype(weights[WEIGHTS_COLUMN]):
        non_numeric_weight_lines = [index + 2 for index in weights[~weights[WEIGHTS_COLUMN].str.isnumeric()].index.tolist()]
        raise InvalidWeightsFile(weights_file, dedent(f"""\
            Found non-numeric weights on the following lines: {non_numeric_weight_lines}
            {WEIGHTS_COLUMN!r} column must be numeric."""))

    if any(weights[WEIGHTS_COLUMN] < 0):
        negative_weight_lines = [index + 2 for index in weights[weights[WEIGHTS_COLUMN] < 0].index.tolist()]
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
