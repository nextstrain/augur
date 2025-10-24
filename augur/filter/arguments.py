"""
Descriptions for augur filter arguments are stored here for reuse by other commands.
"""

from textwrap import dedent
from augur.dates import SUPPORTED_DATE_HELP_TEXT
from augur.filter.io import ACCEPTED_TYPES
from augur.io.metadata import METADATA_DATE_COLUMN
from . import constants


descriptions = {
    "query": dedent("""\
        Filter sequences by attribute. Uses `Pandas DataFrame query syntax`__.
        (e.g., "country == 'Colombia'" or "(country == 'USA' & (division ==
        'Washington'))")

        __ https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query"""),

    "query_columns": dedent(f"""\
        Use alongside query to specify columns and data types in the format
        'column:type', where type is one of
        ({','.join(sorted(ACCEPTED_TYPES))}). Automatic type inference will be
        attempted on all unspecified columns used in the query. Example:
        region:str coverage:float."""),

    "min_date": dedent(f"""\
        Minimal cutoff for date (inclusive). Supported formats:

        """) + SUPPORTED_DATE_HELP_TEXT,

    "max_date": dedent(f"""\
        Maximal cutoff for date (inclusive). Supported formats:

        """) + SUPPORTED_DATE_HELP_TEXT,

    "exclude_ambiguous_dates_by": dedent("""\
        Exclude ambiguous dates by day (e.g., 2020-09-XX), month (e.g.,
        2020-XX-XX), year (e.g., 200X-10-01), or any date fields. An ambiguous
        year makes the corresponding month and day ambiguous, too, even if those
        fields have unambiguous values (e.g., "201X-10-01"). Similarly, an
        ambiguous month makes the corresponding day ambiguous (e.g.,
        "2010-XX-01")."""),

    "exclude": dedent("""\
        File(s) with list of strains to exclude."""),

    "exclude_where": dedent("""\
        Exclude sequences matching these conditions. Ex: "host=rat" or
        "host!=rat". Multiple values are processed as OR (matching any of those
        specified will be excluded), not AND."""),

    "exclude_all": dedent("""\
        Exclude all strains by default. Use this with the include arguments to
        select a specific subset of strains."""),

    "include": dedent("""\
        File(s) with list of strains to include regardless of priorities,
        subsampling, or absence of an entry in sequences."""),

    "include_where": dedent("""\
        Include sequences with these values. ex: host=rat. Multiple values are
        processed as OR (having any of those specified will be included), not
        AND. This rule is applied last and ensures any strains matching these
        rules will be included regardless of priorities, subsampling, or absence
        of an entry in sequences."""),

    "min_length": dedent("""\
        Minimal length of the sequences, only counting standard nucleotide
        characters A, C, G, or T (case-insensitive)."""),

    "max_length": dedent("""\
        Maximum length of the sequences, only counting standard nucleotide
        characters A, C, G, or T (case-insensitive)."""),

    "non_nucleotide": dedent("""\
        Exclude sequences that contain illegal characters."""),

    "group_by": dedent(f"""\
        Grouping columns for subsampling. Notes:

        (1) Grouping by {sorted(constants.GROUP_BY_GENERATED_COLUMNS)} is only
            supported when there is a {METADATA_DATE_COLUMN!r} column in the
            metadata.
        (2) 'week' uses the ISO week numbering system, where a week starts on a
            Monday and ends on a Sunday.
        (3) 'month' and 'week' grouping cannot be used together.
        (4) Custom columns {sorted(constants.GROUP_BY_GENERATED_COLUMNS)} in the
            metadata are ignored for grouping. Please rename them if you want to
            use their values for grouping."""),

    "sequences_per_group": dedent("""\
        Select no more than this number of sequences per category."""),

    "subsample_max_sequences": dedent("""\
        Select no more than this number of sequences (i.e. total sample
        size). Can be used without grouping columns."""),

    "probabilistic_sampling": dedent("""\
        Allow probabilistic sampling during subsampling. This is useful when
        there are more groups than requested sequences. This option only applies
        when a total sample size is provided."""),

    "group_by_weights": dedent("""\
        TSV file defining weights for grouping. Requirements:

        (1) Lines starting with '#' are treated as comment lines.
        (2) The first non-comment line must be a header row.
        (3) There must be a numeric ``weight`` column (weights can take on any
            non-negative values).
        (4) Other columns must be a subset of grouping columns, with
            combinations of values covering all combinations present in the
            metadata.
        (5) This option only applies when grouping columns and a total sample
            size are provided.
        (6) This option can only be used when probabilistic sampling is allowed.

        Notes:

        (1) Any grouping columns absent from this file will be given equal
            weighting across all values *within* groups defined by the other
            weighted columns.
        (2) An entry with the value ``default`` under all columns will be
            treated as the default weight for specific groups present in the
            metadata but missing from the weights file. If there is no default
            weight and the metadata contains rows that are not covered by the
            given weights, augur filter will exit with an error.""")
}
