# Constants set at run time.
RUNTIME_DB_FILE: str = None


# ID column used for all tables defined internally.
ID_COLUMN = 'strain'


# Below are table names, column names, and constant values associated with the database.

# A table representing the original metadata from the user.
METADATA_TABLE = '__augur_filter__metadata'


# A table representing the sequence index, either provided by the user or
# automatically generated.
SEQUENCE_INDEX_TABLE = '__augur_filter__sequence_index'


# A table representing information on where strains are present among multiple
# inputs.
INPUT_SOURCE_TABLE = '__augur_filter__input_source'

# Columns that represent the multiple input sources.
STRAIN_IN_METADATA_COLUMN = '__augur_filter__strain_in_metadata'
STRAIN_IN_SEQUENCES_COLUMN = '__augur_filter__strain_in_sequences'
STRAIN_IN_SEQUENCE_INDEX_COLUMN = '__augur_filter__strain_in_sequence_index'


# A table representing strain priorities for subsampling, either provided by the
# user or automatically generated.
PRIORITIES_TABLE = '__augur_filter__priorities'

# A column for priority scores.
PRIORITY_COLUMN = '__augur_filter__priority'


# A table representing information parsed from the date column in the original
# metadata.
DATE_TABLE = '__augur_filter__metadata_date_info'

# Columns in the date table.
DATE_YEAR_COLUMN = 'year'
DATE_MONTH_COLUMN = 'month'
DATE_DAY_COLUMN = 'day'
DATE_WEEK_COLUMN = 'week'
NUMERIC_DATE_MIN_COLUMN = 'date_min'
NUMERIC_DATE_MAX_COLUMN = 'date_max'

# Generated columns available for --group-by.
# Use sorted() for reproducible output.
GROUP_BY_GENERATED_COLUMNS = {
    DATE_YEAR_COLUMN,
    DATE_MONTH_COLUMN,
    DATE_WEEK_COLUMN,
}


# A table with columns to indicate why and how each row is filtered.
FILTER_REASON_TABLE = '__augur_filter__filter_reason'

# A column for boolean values to indicate excluded strains.
EXCLUDE_COLUMN = '__augur_filter__exclude'

# A column for boolean values to indicate force-included strains.
INCLUDE_COLUMN = '__augur_filter__force_include'

# A column for the filter reason.
FILTER_REASON_COLUMN = 'filter'

# A column for the filter reason's keyword arguments.
FILTER_REASON_KWARGS_COLUMN = 'kwargs'

# A value for the filter reason column to denote exclusion by subsampling.
SUBSAMPLE_FILTER_REASON = 'subsampling'


# A table used only during subsampling that contains the columns necessary for
# grouping.
GROUPING_TABLE = '__augur_filter__grouping'

# A column used when --group-by is not provided to ensure all samples are
# effectively in the same group.
GROUP_BY_DUMMY_COLUMN = '__augur_filter__group_by_placeholder'
GROUP_BY_DUMMY_VALUE = '"dummy"'


# A table used only during subsampling that contains information on group size
# limits based on user-specified parameters.
GROUP_SIZE_LIMITS_TABLE = '__augur_filter__group_size_limits'

# A column for group size limits.
GROUP_SIZE_LIMIT_COLUMN = '__augur_filter__group_size_limit'
