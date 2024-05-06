Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata TSV file for testing.

  $ cat >metadata.tsv <<~~
  > strain	date
  > SEQ_1	2020.0
  > SEQ_2	2020
  > SEQ_3	2020-XX-XX
  > ~~

Test that 2020 is evaluated as 2020-XX-XX.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --min-date 2020-02-01 \
  >  --output-strains filtered_strains.txt 2>/dev/null
  $ sort filtered_strains.txt
  SEQ_2
  SEQ_3

Test that 2020.0, 2020, and 2020-XX-XX all pass --min-date 2019

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --min-date 2019 \
  >  --output-strains filtered_strains.txt 2>/dev/null
  $ sort filtered_strains.txt
  SEQ_1
  SEQ_2
  SEQ_3
