Setup

  $ source "$TESTDIR"/_setup.sh

SEQ1 and SEQ2 translate to week=(2003, 1).
SEQ3 and SEQ4 translate to week=(2004, 1).
These should be in separate groups.

  $ cat >metadata.tsv <<~~
  > strain	date
  > SEQ1	2003-01-01
  > SEQ2	2003-01-02
  > SEQ3	2003-12-30
  > SEQ4	2003-12-31
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by week \
  >   --sequences-per-group 1 \
  >   --subsample-seed 0 \
  >   --output-metadata metadata-filtered.tsv
  2 strains were dropped during filtering
  	2 were dropped because of subsampling criteria
  2 strains passed all filters
  $ cat metadata-filtered.tsv
  strain	date
  SEQ1	2003-01-01
  SEQ3	2003-12-30

ISO year from 'week' takes precedence over 'year'.

  $ cat >metadata.tsv <<~~
  > strain	date
  > SEQ1	2003-12-30
  > SEQ2	2003-12-31
  > SEQ3	2004-01-01
  > SEQ4	2004-01-02
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by year week \
  >   --sequences-per-group 1 \
  >   --subsample-seed 0 \
  >   --output-metadata metadata-filtered.tsv
  WARNING: 'year' grouping will be ignored since 'week' includes ISO year.
  3 strains were dropped during filtering
  	3 were dropped because of subsampling criteria
  1 strain passed all filters
  $ cat metadata-filtered.tsv
  strain	date
  SEQ1	2003-12-30
