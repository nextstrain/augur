Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata file for testing.

  $ cat >metadata.tsv <<~~
  > strain	region	quality
  > SEQ_1	Europe	good
  > SEQ_2	Europe	bad
  > SEQ_3	Asia	good
  > ~~

Multiple --query values are combined: only sequences matching every query are
kept.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "region == 'Europe'" "quality == 'good'" \
  >  --output-strains filtered_strains.txt 2>/dev/null

  $ cat filtered_strains.txt
  SEQ_1

This is equivalent to a single query joining conditions with `&`.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "(region == 'Europe') & (quality == 'good')" \
  >  --output-strains filtered_strains.txt 2>/dev/null

  $ cat filtered_strains.txt
  SEQ_1

Queries are shown separately in the final report.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "region == 'Europe'" "quality == 'good'" \
  >  --output-strains filtered_strains.txt 2>&1 >/dev/null
  2 strains were dropped during filtering
  \t1 was filtered out by the query: "region == 'Europe'" (esc)
  \t1 was filtered out by the query: "quality == 'good'" (esc)
  1 strain passed all filters
