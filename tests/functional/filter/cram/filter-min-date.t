Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata TSV file for testing.

  $ cat >metadata.tsv <<~~
  > strain	date
  > SEQ_1	2020-02-XX
  > SEQ_2	2020-02-26
  > SEQ_3	2020-02-25
  > ~~

Test that --min-date is inclusive.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --min-date 2020-02-26 \
  >  --output-strains filtered_strains.txt 2>/dev/null
  $ sort filtered_strains.txt
  SEQ_1
  SEQ_2
