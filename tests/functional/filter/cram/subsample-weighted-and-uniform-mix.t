Setup

  $ source "$TESTDIR"/_setup.sh

Generate metadata file with 250 rows.

  $ echo "strain	date	location" > metadata.tsv
  $ for i in $(seq 1 50); do
  >     echo "2000A_$i	2000	A" >> metadata.tsv
  >     echo "2000B_$i	2000	B" >> metadata.tsv
  >     echo "2001A_$i	2001	A" >> metadata.tsv
  >     echo "2001B_$i	2001	B" >> metadata.tsv
  >     echo "2002B_$i	2002	B" >> metadata.tsv
  > done

Weight locations A:B as 2:1. This is reflected in target_group_sizes.tsv below.

  $ cat >weights.tsv <<~~
  > location	weight
  > A	2
  > B	1
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by location \
  >   --group-by-weights weights.tsv \
  >   --subsample-max-sequences 100 \
  >   --subsample-seed 0 \
  >   --output-group-by-sizes target_group_sizes.tsv \
  >   --output-strains strains.txt 2>/dev/null

  $ cat target_group_sizes.tsv
  location	target_size
  A	67
  B	33

Using 1:1 weights is similarly straightforward, with 50 sequences from each location.

  $ cat >weights.tsv <<~~
  > location	weight
  > A	1
  > B	1
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by location \
  >   --group-by-weights weights.tsv \
  >   --subsample-max-sequences 100 \
  >   --subsample-seed 0 \
  >   --output-group-by-sizes target_group_sizes.tsv \
  >   --output-strains strains.txt 2>/dev/null

  $ cat target_group_sizes.tsv
  location	target_size
  A	50
  B	50

Keep the 1:1 location weighting, but add uniform sampling on year.
Since there are no available sequences in the group (2002,A), we end up with
more sequences from location B relative to A.

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by year location \
  >   --group-by-weights weights.tsv \
  >   --subsample-max-sequences 100 \
  >   --subsample-seed 0 \
  >   --output-group-by-sizes target_group_sizes.tsv \
  >   --output-strains strains.txt 2>/dev/null

  $ cat target_group_sizes.tsv
  year	location	target_size
  2000	A	20
  2000	B	20
  2001	A	20
  2001	B	20
  2002	B	20

If a single sequence is added for group (2002,A), the weighting now appears
"equal" among all years and locations.

However, there is only 1 sequence available in (2002,A), much lower than the
requested 17, so the total number of sequences outputted is lower than requested.

  $ echo "2002A_1	2002	A" >> metadata.tsv

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by year location \
  >   --group-by-weights weights.tsv \
  >   --subsample-max-sequences 100 \
  >   --subsample-seed 0 \
  >   --output-group-by-sizes target_group_sizes.tsv \
  >   --output-strains strains.txt 2>/dev/null

  $ cat target_group_sizes.tsv
  year	location	target_size
  2000	A	17
  2000	B	16
  2001	A	16
  2001	B	16
  2002	A	17
  2002	B	17

  $ wc -l strains.txt
  \s*83 .* (re)
