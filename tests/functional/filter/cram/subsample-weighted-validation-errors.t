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

  $ cat >weights.tsv <<~~
  > location	weight
  > A	2
  > B	1
  > ~~

When --group-by-weights is specified, all columns must be provided in
--group-by.

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

--output-group-by-sizes is only available for --group-by-weights.

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by location \
  >   --subsample-max-sequences 6 \
  >   --subsample-seed 0 \
  >   --output-group-by-sizes sizes.tsv \
  >   --output-strains strains.txt
  ERROR: --output-group-by-sizes is only available for --group-by-weights. It may be added to other sampling methods in the future - see <https://github.com/nextstrain/augur/issues/1590>
  [2]

--group-by-weights cannot be used with --no-probabilistic-sampling.

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by location \
  >   --group-by-weights weights.tsv \
  >   --no-probabilistic-sampling \
  >   --subsample-max-sequences 6 \
  >   --subsample-seed 0 \
  >   --output-strains strains.txt
  ERROR: --group-by-weights cannot be used with --no-probabilistic-sampling.
  [2]
