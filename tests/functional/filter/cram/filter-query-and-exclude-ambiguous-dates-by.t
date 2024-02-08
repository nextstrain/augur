Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata TSV file for testing.

  $ cat >metadata.tsv <<~~
  > strain	date	region
  > SEQ_1	2020	Asia
  > SEQ_2	2020	Asia
  > SEQ_3	2020	Asia
  > SEQ_4	2020	North America
  > ~~

Confirm that `--exclude-ambiguous-dates-by` works for all year only ambiguous dates.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query 'region=="Asia"' \
  >  --exclude-ambiguous-dates-by any \
  >  --empty-output-reporting silent \
  >  --output-strains filtered_strains.txt
  4 strains were dropped during filtering
  	1 was filtered out by the query: "region=="Asia""
  	3 were dropped because of their ambiguous date in any
  0 strains passed all filters
