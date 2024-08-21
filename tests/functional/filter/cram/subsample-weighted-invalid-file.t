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

Weights must be non-negative.

  $ cat >weights.tsv <<~~
  > location	weight
  > A	2
  > B	1
  > C	-1
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by location \
  >   --group-by-weights weights.tsv \
  >   --subsample-max-sequences 6 \
  >   --subsample-seed 0 \
  >   --output-strains strains.txt
  Sampling with weights defined by weights.tsv.
  ERROR: Bad weights file 'weights.tsv'.
  Found negative weights on the following lines: [4]
  'weight' column must be non-negative.
  [2]

Weights must be numeric.

  $ cat >weights.tsv <<~~
  > location	weight
  > A	yes
  > B	1
  > C	no
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by location \
  >   --group-by-weights weights.tsv \
  >   --subsample-max-sequences 6 \
  >   --subsample-seed 0 \
  >   --output-strains strains.txt
  Sampling with weights defined by weights.tsv.
  ERROR: Bad weights file 'weights.tsv'.
  Found non-numeric weights on the following lines: [2, 4]
  'weight' column must be numeric.
  [2]

Weights file cannot be empty.

  $ cat >weights.tsv <<~~
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by location \
  >   --group-by-weights weights.tsv \
  >   --subsample-max-sequences 6 \
  >   --subsample-seed 0 \
  >   --output-strains strains.txt
  ERROR: Bad weights file 'weights.tsv'.
  File is empty.
  [2]
