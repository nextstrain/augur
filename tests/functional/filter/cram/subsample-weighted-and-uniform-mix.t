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

  $ cat >weights-A2B1.tsv <<~~
  > location	weight
  > A	2
  > B	1
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by location \
  >   --group-by-weights weights-A2B1.tsv \
  >   --subsample-max-sequences 100 \
  >   --subsample-seed 0 \
  >   --output-group-by-sizes target_group_sizes.tsv \
  >   --output-metadata filtered.tsv 2>/dev/null

  $ cat target_group_sizes.tsv
  location	weight	_augur_filter_target_size	_augur_filter_input_size	_augur_filter_output_size
  A	2	67	100	67
  B	1	33	150	33

There are also enough rows per group that the output metadata directly reflects
the target group sizes.

  $ cat filtered.tsv | tail -n +2 | cut -f3 | sort | uniq -c
  \s*67 A (re)
  \s*33 B (re)

Using 1:1 weights is similarly straightforward, with 50 sequences from each location.

  $ cat >weights-A1B1.tsv <<~~
  > location	weight
  > A	1
  > B	1
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by location \
  >   --group-by-weights weights-A1B1.tsv \
  >   --subsample-max-sequences 100 \
  >   --subsample-seed 0 \
  >   --output-group-by-sizes target_group_sizes.tsv \
  >   --output-strains strains.txt 2>/dev/null

  $ cat target_group_sizes.tsv
  location	weight	_augur_filter_target_size	_augur_filter_input_size	_augur_filter_output_size
  A	1	50	100	50
  B	1	50	150	50

Keep the 1:1 location weighting, but add uniform sampling on year.
The uniform sampling happens "within" each weighted column value, so the 1:1
location weighting is reflected even though there is an imbalance in years
available per location.

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by year location \
  >   --group-by-weights weights-A1B1.tsv \
  >   --subsample-max-sequences 100 \
  >   --subsample-seed 0 \
  >   --output-group-by-sizes target_group_sizes.tsv \
  >   --output-strains strains.txt 2>/dev/null

  $ cat target_group_sizes.tsv
  year	location	weight	_augur_filter_target_size	_augur_filter_input_size	_augur_filter_output_size
  2000	A	0.5	25	50	25
  2000	B	0.3333333333333333	16	50	16
  2001	A	0.5	25	50	25
  2001	B	0.3333333333333333	16	50	16
  2002	B	0.3333333333333333	17	50	17

If a single sequence is added for group (2002,A), the weighting now appears
"equal" among all years and locations.

However, there is only 1 sequence available in (2002,A), much lower than the
requested 17, so the total number of sequences outputted is lower than requested.

  $ echo "2002A_1	2002	A" >> metadata.tsv

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by year location \
  >   --group-by-weights weights-A1B1.tsv \
  >   --subsample-max-sequences 100 \
  >   --subsample-seed 0 \
  >   --output-group-by-sizes target_group_sizes.tsv \
  >   --output-strains strains.txt 2>/dev/null

  $ cat target_group_sizes.tsv
  year	location	weight	_augur_filter_target_size	_augur_filter_input_size	_augur_filter_output_size
  2000	A	0.3333333333333333	17	50	17
  2000	B	0.3333333333333333	16	50	16
  2001	A	0.3333333333333333	16	50	16
  2001	B	0.3333333333333333	16	50	16
  2002	A	0.3333333333333333	17	1	1
  2002	B	0.3333333333333333	17	50	17

  $ wc -l strains.txt
  \s*83 .* (re)
