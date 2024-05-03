Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata file for testing.

  $ cat >metadata.tsv <<~~
  > strain	column
  > SEQ_1	value1
  > SEQ_2	value2
  > SEQ_3	value3
  > ~~

'column' should be query-able using the `.str` accessor.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "column.str.startswith('value')" \
  >  --output-strains filtered_strains.txt 2>/dev/null

  $ sort filtered_strains.txt
  SEQ_1
  SEQ_2
  SEQ_3
