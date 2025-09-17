Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata TSV file for testing.

  $ cat >metadata.tsv <<~~
  > strain	date
  > SEQ_1	2019-XX-XX
  > SEQ_2	2019-12-31
  > SEQ_3	2020-01-01
  > ~~

Test that --max-date is inclusive even with ambiguity.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --max-date 2019 \
  >  --output-strains filtered_strains.txt 2>/dev/null
  $ sort filtered_strains.txt
  SEQ_1
  SEQ_2
