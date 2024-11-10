Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata TSV file for testing.

  $ cat >metadata.tsv <<~~
  > strain	date
  > SEQ_1	2020
  > SEQ_2	2020
  > SEQ_3	2020
  > ~~

Confirm that `--exclude-ambiguous-dates-by` works for all year only ambiguous dates.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --exclude-ambiguous-dates-by any \
  >  --empty-output-reporting silent \
  >  --output-strains filtered_strains.txt
  3 strains were dropped during filtering
  	3 were dropped because of their ambiguous date in any
  0 strains passed all filters
