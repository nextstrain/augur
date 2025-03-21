Setup

  $ source "$TESTDIR"/_setup.sh

Set up files.

  $ cat >metadata.tsv <<~~
  > strain	date
  > SEQ1	2001-01-01
  > SEQ2	2001-01-01
  > SEQ3	2002-01-01
  > SEQ4	2002-01-01
  > SEQ5	2003-01-01
  > SEQ6	2003-01-01
  > ~~

  $ cat >weights.tsv <<~~
  > year	weight
  > 2001	2
  > 2002	2
  > default	1
  > ~~

Subsample with year weights.

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by year \
  >   --group-by-weights weights.tsv \
  >   --subsample-max-sequences 5 \
  >   --subsample-seed 0 \
  >   --output-strains strains.txt
  Sampling with weights defined by weights.tsv.
  WARNING: The input metadata contains these values under the following columns that are not directly covered by 'weights.tsv':
  - 'year': ['2003']
  The default weight of 1 will be used for all groups defined by those values.
  NOTE: Skipping 1 group due to lack of entries in metadata.
  1 strain was dropped during filtering
  	1 was dropped because of subsampling criteria
  5 strains passed all filters
