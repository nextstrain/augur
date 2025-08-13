"""
Descriptions for augur filter arguments are stored here for reuse by other commands.
"""

from augur.dates import SUPPORTED_DATE_HELP_TEXT
from augur.filter.io import ACCEPTED_TYPES
from augur.io.metadata import METADATA_DATE_COLUMN
from . import constants


descriptions = {
    "query": """Filter samples by attribute.
        Uses Pandas Dataframe querying, see https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query for syntax.
        (e.g., --query "country == 'Colombia'" or --query "(country == 'USA' & (division == 'Washington'))")""",

    "query_columns": f"""
        Use alongside --query to specify columns and data types in the format 'column:type', where type is one of ({','.join(sorted(ACCEPTED_TYPES))}).
        Automatic type inference will be attempted on all unspecified columns used in the query.
        Example: region:str coverage:float.
    """,

    "min_date": f"minimal cutoff for date, the cutoff date is inclusive; may be specified as: {SUPPORTED_DATE_HELP_TEXT}",

    "max_date": f"maximal cutoff for date, the cutoff date is inclusive; may be specified as: {SUPPORTED_DATE_HELP_TEXT}",

    "exclude_ambiguous_dates_by": 'Exclude ambiguous dates by day (e.g., 2020-09-XX), month (e.g., 2020-XX-XX), year (e.g., 200X-10-01), or any date fields. An ambiguous year makes the corresponding month and day ambiguous, too, even if those fields have unambiguous values (e.g., "201X-10-01"). Similarly, an ambiguous month makes the corresponding day ambiguous (e.g., "2010-XX-01").',

    "exclude": "file(s) with list of strains to exclude",

    "exclude_where": "Exclude samples matching these conditions. Ex: \"host=rat\" or \"host!=rat\". Multiple values are processed as OR (matching any of those specified will be excluded), not AND",

    "exclude_all": "exclude all strains by default. Use this with the include arguments to select a specific subset of strains.",

    "include": "file(s) with list of strains to include regardless of priorities, subsampling, or absence of an entry in --sequences.",

    "include_where": """
        Include samples with these values. ex: host=rat. Multiple values are
        processed as OR (having any of those specified will be included), not
        AND. This rule is applied last and ensures any strains matching these
        rules will be included regardless of priorities, subsampling, or absence
        of an entry in --sequences.""",

    "min_length": "minimal length of the sequences, only counting standard nucleotide characters A, C, G, or T (case-insensitive)",

    "max_length": "maximum length of the sequences, only counting standard nucleotide characters A, C, G, or T (case-insensitive)",

    "non_nucleotide": "exclude sequences that contain illegal characters",

    "group_by": f"""
        categories with respect to subsample.
        Notes:
        (1) Grouping by {sorted(constants.GROUP_BY_GENERATED_COLUMNS)} is only supported when there is a {METADATA_DATE_COLUMN!r} column in the metadata.
        (2) 'week' uses the ISO week numbering system, where a week starts on a Monday and ends on a Sunday.
        (3) 'month' and 'week' grouping cannot be used together.
        (4) Custom columns {sorted(constants.GROUP_BY_GENERATED_COLUMNS)} in the metadata are ignored for grouping. Please rename them if you want to use their values for grouping.""",

    "sequences_per_group": "subsample to no more than this number of sequences per category",

    "subsample_max_sequences": "subsample to no more than this number of sequences; can be used without the group_by argument",

    "probabilistic_sampling": "Allow probabilistic sampling during subsampling. This is useful when there are more groups than requested sequences. This option only applies when `--subsample-max-sequences` is provided.",

    "group_by_weights": """
        TSV file defining weights for grouping. Requirements:

        (1) Lines starting with '#' are treated as comment lines.
        (2) The first non-comment line must be a header row.
        (3) There must be a numeric ``weight`` column (weights can take on any
            non-negative values).
        (4) Other columns must be a subset of columns used in ``--group-by``,
            with combinations of values covering all combinations present in the
            metadata.
        (5) This option only applies when ``--group-by`` and
            ``--subsample-max-sequences`` are provided.
        (6) This option cannot be used with ``--no-probabilistic-sampling``.

        Notes:

        (1) Any ``--group-by`` columns absent from this file will be given equal
            weighting across all values *within* groups defined by the other
            weighted columns.
        (2) An entry with the value ``default`` under all columns will be
            treated as the default weight for specific groups present in the
            metadata but missing from the weights file. If there is no default
            weight and the metadata contains rows that are not covered by the
            given weights, augur filter will exit with an error.
    """
}
