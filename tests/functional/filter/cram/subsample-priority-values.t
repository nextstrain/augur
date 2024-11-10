Setup

  $ source "$TESTDIR"/_setup.sh

Create files for testing.

  $ cat >metadata.tsv <<~~
  > strain	quality
  > SEQ_1	good
  > SEQ_2	good
  > SEQ_3	good
  > SEQ_4	good
  > SEQ_5	good
  > ~~

  $ cat >priorities.tsv <<~~
  > SEQ_1	5
  > SEQ_2	6
  > SEQ_3	8
  > SEQ_4	-100
  > ~~

Subsample 5 strains. Note that SEQ_5 should not be dropped even though it does
not have a priority score.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --priority priorities.tsv \
  >  --subsample-max-sequences 5 \
  >  --output-strains filtered_strains.txt 2>/dev/null
  $ sort filtered_strains.txt
  SEQ_1
  SEQ_2
  SEQ_3
  SEQ_4
  SEQ_5

Subsample 1 less strain. SEQ_5 should be dropped since it does not have a
priority score, which is effectively lowest priority (even compared to negative
values).

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --priority priorities.tsv \
  >  --subsample-max-sequences 4 \
  >  --output-strains filtered_strains.txt 2>/dev/null
  $ sort filtered_strains.txt
  SEQ_1
  SEQ_2
  SEQ_3
  SEQ_4

Subsample 1 less strain. SEQ_4 should now be dropped.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --priority priorities.tsv \
  >  --subsample-max-sequences 3 \
  >  --output-strains filtered_strains.txt 2>/dev/null
  $ sort filtered_strains.txt
  SEQ_1
  SEQ_2
  SEQ_3

Subsample 1 less strain. SEQ_1 should now be dropped.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --priority priorities.tsv \
  >  --subsample-max-sequences 2 \
  >  --output-strains filtered_strains.txt 2>/dev/null
  $ sort filtered_strains.txt
  SEQ_2
  SEQ_3

Subsample 1 less strain. SEQ_1 should now be dropped.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --priority priorities.tsv \
  >  --subsample-max-sequences 1 \
  >  --output-strains filtered_strains.txt 2>/dev/null
  $ sort filtered_strains.txt
  SEQ_3
