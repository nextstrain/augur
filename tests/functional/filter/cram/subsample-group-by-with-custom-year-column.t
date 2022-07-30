Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Create a metadata file with a custom year column

  $ cat >$TMP/metadata-year-column.tsv <<~~
  > strain	date	year	month
  > SEQ1	2021-01-01	odd	January
  > SEQ2	2021-01-02	odd	January
  > SEQ3	2022-01-01	even	January
  > SEQ4	2022-01-02	even	January
  > SEQ5	2022-02-02	even	February
  > ~~

Pandas engine
-------------

Group by generated year column, and ensure all original columns are still in the final output.

  $ ${AUGUR} filter \
  >  --metadata $TMP/metadata-year-column.tsv \
  >  --group-by year \
  >  --sequences-per-group 1 \
  >  --subsample-seed 0 \
  >  --output-metadata "$TMP/filtered_metadata.tsv" > /dev/null
  WARNING: `--group-by year` uses the generated year value from the 'date' column. The custom 'year' column in the metadata is ignored for grouping purposes.
  $ cat "$TMP/filtered_metadata.tsv"
  strain\tdate\tyear\tmonth (esc)
  SEQ1\t2021-01-01\todd\tJanuary (esc)
  SEQ5\t2022-02-02\teven\tFebruary (esc)

Group by generated year and month columns, and ensure all original columns are still in the final output.

  $ ${AUGUR} filter \
  >  --metadata $TMP/metadata-year-column.tsv \
  >  --group-by year month \
  >  --sequences-per-group 1 \
  >  --subsample-seed 0 \
  >  --output-metadata "$TMP/filtered_metadata.tsv" > /dev/null
  WARNING: `--group-by year` uses the generated year value from the 'date' column. The custom 'year' column in the metadata is ignored for grouping purposes.
  WARNING: `--group-by month` uses the generated month value from the 'date' column. The custom 'month' column in the metadata is ignored for grouping purposes.
  $ cat "$TMP/filtered_metadata.tsv"
  strain\tdate\tyear\tmonth (esc)
  SEQ1\t2021-01-01\todd\tJanuary (esc)
  SEQ3\t2022-01-01\teven\tJanuary (esc)
  SEQ5\t2022-02-02\teven\tFebruary (esc)

SQLite engine
-------------

Group by generated year column, and ensure all original columns are still in the final output.

  $ ${AUGUR} filter --engine sqlite \
  >  --metadata $TMP/metadata-year-column.tsv \
  >  --group-by year \
  >  --sequences-per-group 1 \
  >  --subsample-seed 0 \
  >  --output-metadata "$TMP/filtered_metadata.tsv" > /dev/null
  WARNING: `--group-by year` uses the generated year value from the 'date' column. The custom 'year' column in the metadata is ignored for grouping purposes.
  $ cat "$TMP/filtered_metadata.tsv"
  strain\tdate\tyear\tmonth (esc)
  SEQ1\t2021-01-01\todd\tJanuary (esc)
  SEQ5\t2022-02-02\teven\tFebruary (esc)

Group by generated year and month columns, and ensure all original columns are still in the final output.

  $ ${AUGUR} filter --engine sqlite \
  >  --metadata $TMP/metadata-year-column.tsv \
  >  --group-by year month \
  >  --sequences-per-group 1 \
  >  --subsample-seed 0 \
  >  --output-metadata "$TMP/filtered_metadata.tsv" > /dev/null
  WARNING: `--group-by year` uses the generated year value from the 'date' column. The custom 'year' column in the metadata is ignored for grouping purposes.
  WARNING: `--group-by month` uses the generated month value from the 'date' column. The custom 'month' column in the metadata is ignored for grouping purposes.
  $ cat "$TMP/filtered_metadata.tsv"
  strain\tdate\tyear\tmonth (esc)
  SEQ1\t2021-01-01\todd\tJanuary (esc)
  SEQ3\t2022-01-01\teven\tJanuary (esc)
  SEQ5\t2022-02-02\teven\tFebruary (esc)
