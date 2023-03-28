# Shared variables set at run time.
# TODO: Remove these with the database implementation.
sequence_index = None
sequence_strains = None
metadata_strains = None
valid_strains = None
all_sequences_to_include = None
filter_counts = None
num_excluded_subsamp = None

# Generated date columns.
DATE_YEAR_COLUMN = 'year'
DATE_MONTH_COLUMN = 'month'
DATE_WEEK_COLUMN = 'week'

# Generated columns available for --group-by.
# Use sorted() for reproducible output.
GROUP_BY_GENERATED_COLUMNS = {
    DATE_YEAR_COLUMN,
    DATE_MONTH_COLUMN,
    DATE_WEEK_COLUMN,
}
