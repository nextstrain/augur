Setup

  $ source "$TESTDIR"/_setup.sh

Create a metadata file with a custom year column

  $ cat >metadata-year-column.tsv <<~~
  > strain	date	year	month
  > SEQ1	2021-01-01	odd	January
  > SEQ2	2021-01-02	odd	January
  > SEQ3	2022-01-01	even	January
  > SEQ4	2022-01-02	even	January
  > SEQ5	2022-02-02	even	February
  > ~~

Group by generated year column, and ensure all original columns are still in the final output.

  $ ${AUGUR} filter \
  >  --metadata metadata-year-column.tsv \
  >  --group-by year \
  >  --sequences-per-group 1 \
  >  --subsample-seed 0 \
  >  --output-metadata filtered_metadata.tsv > /dev/null
  WARNING: `--group-by year` uses a generated year value from the 'date' column. The custom 'year' column in the metadata is ignored for grouping purposes.
  3 strains were dropped during filtering
  	3 were dropped because of subsampling criteria
  2 strains passed all filters
  $ cat filtered_metadata.tsv
  strain\tdate\tyear\tmonth (esc)
  SEQ1\t2021-01-01\todd\tJanuary (esc)
  SEQ5\t2022-02-02\teven\tFebruary (esc)

Group by generated year and month columns, and ensure all original columns are still in the final output.

  $ ${AUGUR} filter \
  >  --metadata metadata-year-column.tsv \
  >  --group-by year month \
  >  --sequences-per-group 1 \
  >  --subsample-seed 0 \
  >  --output-metadata filtered_metadata.tsv > /dev/null
  WARNING: `--group-by month` uses a generated month value from the 'date' column. The custom 'month' column in the metadata is ignored for grouping purposes.
  WARNING: `--group-by year` uses a generated year value from the 'date' column. The custom 'year' column in the metadata is ignored for grouping purposes.
  2 strains were dropped during filtering
  	2 were dropped because of subsampling criteria
  3 strains passed all filters
  $ cat filtered_metadata.tsv
  strain\tdate\tyear\tmonth (esc)
  SEQ1\t2021-01-01\todd\tJanuary (esc)
  SEQ3\t2022-01-01\teven\tJanuary (esc)
  SEQ5\t2022-02-02\teven\tFebruary (esc)
