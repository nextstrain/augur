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

Sampling with location weights only.

  $ cat >weights.tsv <<~~
  > location	weight
  > A	2
  > B	1
  > C	3
  > ~~

This should take 4 from location A and 2 from location B. The weight for
location C is ignored because there are no corresponding rows in the metadata.

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by location \
  >   --group-by-weights weights.tsv \
  >   --subsample-max-sequences 6 \
  >   --subsample-seed 0 \
  >   --output-strains strains.txt
  Sampling with weights defined by weights.tsv.
  Skipping 1 group due to lack of entries in metadata.
  2 strains were dropped during filtering
  	2 were dropped because of subsampling criteria
  6 strains passed all filters

  $ cat strains.txt
  SEQ1
  SEQ2
  SEQ5
  SEQ6
  SEQ7
  SEQ8

Sampling with weights on location and uniform sampling on date (--group-by
month) should work.

  $ cat >weights.tsv <<~~
  > location	weight
  > A	2
  > B	1
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by location month \
  >   --group-by-weights weights.tsv \
  >   --subsample-max-sequences 6 \
  >   --subsample-seed 0 \
  >   --output-strains strains.txt
  Sampling with weights defined by weights.tsv.
  Skipping 2 groups due to lack of entries in metadata.
  2 strains were dropped during filtering
  	2 were dropped because of subsampling criteria
  6 strains passed all filters

Sampling with incomplete weights should raise an error.

  $ cat >weights.tsv <<~~
  > location	weight
  > A	2
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by location \
  >   --group-by-weights weights.tsv \
  >   --subsample-max-sequences 6 \
  >   --subsample-seed 0 \
  >   --output-strains strains.txt
  Sampling with weights defined by weights.tsv.
  ERROR: The following groups appear in the metadata but are missing from the weights file:
  	location='B'
  [2]

When --group-by-weights is specified, all columns must be provided in
--group-by.

  $ cat >weights.tsv <<~~
  > location	weight
  > A	2
  > B	1
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by-weights weights.tsv \
  >   --subsample-max-sequences 6 \
  >   --subsample-seed 0 \
  >   --output-strains strains.txt
  ERROR: Columns in --group-by-weights must be a subset of columns provided in --group-by.
  [2]

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by month \
  >   --group-by-weights weights.tsv \
  >   --subsample-max-sequences 6 \
  >   --subsample-seed 0 \
  >   --output-strains strains.txt
  ERROR: Columns in --group-by-weights must be a subset of columns provided in --group-by.
  [2]

Negative weights are not allowed.

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
  ERROR: Weights must be non-negative, but found negative weights:   location  weight
  2        C      -1
  [2]
