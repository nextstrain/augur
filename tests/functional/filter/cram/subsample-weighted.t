Setup

  $ source "$TESTDIR"/_setup.sh

Set up files.

  $ cat >metadata.tsv <<~~
  > strain	custom_month	location
  > SEQ1	2000-01	A
  > SEQ2	2000-01	A
  > SEQ3	2000-01	B
  > SEQ4	2000-01	B
  > ~~

  $ cat >weights.tsv <<~~
  > custom_month	location	weight
  > 2000-01	A	1
  > 2000-01	B	2
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --weights weights.tsv \
  >   --subsample-max-sequences 3 \
  >   --subsample-seed 0 \
  >   --output-strains strains.txt
