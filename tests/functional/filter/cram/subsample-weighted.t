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
  NOTE: Skipping 1 group due to lack of entries in metadata.
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
  NOTE: Skipping 2 groups due to lack of entries in metadata.
  NOTE: Weights were not provided for the column 'month'. Using equal weights across values in that column.
  2 strains were dropped during filtering
  	2 were dropped because of subsampling criteria
  6 strains passed all filters

Sampling with incomplete weights should show an error.

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
  ERROR: The input metadata contains these values under the following columns that are not covered by 'weights.tsv':
  - 'location': ['B']
  To fix this, either:
  (1) specify weights explicitly - add entries to 'weights.tsv' for the values above, or
  (2) specify a default weight - add an entry to 'weights.tsv' with the value 'default' for all columns
  [2]

Re-running with a default weight shows a warning and continues.

  $ cat >weights.tsv <<~~
  > location	weight
  > A	2
  > default	1
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by month location \
  >   --group-by-weights weights.tsv \
  >   --subsample-max-sequences 6 \
  >   --subsample-seed 0 \
  >   --output-strains strains.txt
  Sampling with weights defined by weights.tsv.
  WARNING: The input metadata contains these values under the following columns that are not directly covered by 'weights.tsv':
  - 'location': ['B']
  The default weight of 1 will be used for all groups defined by those values.
  NOTE: Skipping 4 groups due to lack of entries in metadata.
  NOTE: Weights were not provided for the column 'month'. Using equal weights across values in that column.
  2 strains were dropped during filtering
  	2 were dropped because of subsampling criteria
  6 strains passed all filters

To specify a default weight, the value 'default' must be set for all weighted columns.

  $ cat >weights.tsv <<~~
  > location	month	weight
  > A	2000-01	2
  > A	2000-02	2
  > default		1
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by month location \
  >   --group-by-weights weights.tsv \
  >   --subsample-max-sequences 6 \
  >   --subsample-seed 0 \
  >   --output-strains strains.txt
  Sampling with weights defined by weights.tsv.
  ERROR: The input metadata contains these values under the following columns that are not covered by 'weights.tsv':
  - 'location': ['B']
  - 'month': ['2000-01', '2000-03']
  To fix this, either:
  (1) specify weights explicitly - add entries to 'weights.tsv' for the values above, or
  (2) specify a default weight - add an entry to 'weights.tsv' with the value 'default' for all columns
  [2]

  $ cat >weights.tsv <<~~
  > location	month	weight
  > A	2000-01	2
  > A	2000-02	2
  > default	default	1
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by month location \
  >   --group-by-weights weights.tsv \
  >   --subsample-max-sequences 6 \
  >   --subsample-seed 0 \
  >   --output-strains strains.txt
  Sampling with weights defined by weights.tsv.
  WARNING: The input metadata contains these values under the following columns that are not directly covered by 'weights.tsv':
  - 'location': ['B']
  - 'month': ['2000-01', '2000-03']
  The default weight of 1 will be used for all groups defined by those values.
  NOTE: Skipping 1 group due to lack of entries in metadata.
  2 strains were dropped during filtering
  	2 were dropped because of subsampling criteria
  6 strains passed all filters
