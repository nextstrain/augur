Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata TSV file for testing.

  $ cat >metadata.tsv <<~~
  > strain	location	quality
  > SEQ_1	colorado	good
  > SEQ_2	colorado	bad
  > SEQ_3	nevada	good
  > ~~
  $ cat >include.txt <<~~
  > SEQ_3
  > ~~

Test that --include_where still works with filtering on query.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "quality=='good' & location=='colorado'" \
  >  --include include.txt \
  >  --output-strains filtered_strains.txt 2>/dev/null
  $ sort filtered_strains.txt
  SEQ_1
  SEQ_3
