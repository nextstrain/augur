Setup

  $ source "$TESTDIR"/_setup.sh

Set up files.

  $ cat >metadata.tsv <<~~
  > strain	date	location
  > SEQ1	2000-01-01	A
  > SEQ2	2000-01-02	A
  > SEQ3	2000-01-03	B
  > SEQ4	2000-01-04	B
  > SEQ5	2000-02-01	A
  > SEQ6	2000-02-02	A
  > SEQ7	2000-03-01	B
  > SEQ8	2000-03-02	B
  > ~~

Comments in the weights file are valid.

  $ cat >weights.tsv <<~~
  > # This is a comment
  > ## So is this
  > location	weight
  > A	2
  > B	1
  > C	3
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by location \
  >   --group-by-weights weights.tsv \
  >   --subsample-max-sequences 6 \
  >   --output-strains strains.txt
  Sampling with weights defined by weights.tsv.
  NOTE: Skipping 1 group due to lack of entries in metadata.
  2 strains were dropped during filtering
  	2 were dropped because of subsampling criteria
  6 strains passed all filters
