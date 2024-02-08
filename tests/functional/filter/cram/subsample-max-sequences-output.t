Setup

  $ source "$TESTDIR"/_setup.sh

Create a metadata file with 2 groups with sizes (1,5).

  $ cat >metadata.tsv <<~~
  > strain	group
  > SEQ1-A	A
  > SEQ2-B	B
  > SEQ3-B	B
  > SEQ4-B	B
  > SEQ5-B	B
  > SEQ6-B	B
  > ~~

Request a maximum of 9 strains. All 6 strains should be used, even though this means groups are sampled unevenly.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --group-by group \
  >  --subsample-max-sequences 9 \
  >  --output-strains filtered_strains.txt
  Sampling at 9 per group.
  0 strains were dropped during filtering
  	0 were dropped because of subsampling criteria
  6 strains passed all filters
  $ sort filtered_strains.txt
  SEQ1-A
  SEQ2-B
  SEQ3-B
  SEQ4-B
  SEQ5-B
  SEQ6-B

Request a maximum of 6 strains. All 6 strains should be used, even though this means groups are sampled unevenly.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --group-by group \
  >  --subsample-max-sequences 6 \
  >  --output-strains filtered_strains.txt
  Sampling at 6 per group.
  0 strains were dropped during filtering
  	0 were dropped because of subsampling criteria
  6 strains passed all filters
  $ sort filtered_strains.txt
  SEQ1-A
  SEQ2-B
  SEQ3-B
  SEQ4-B
  SEQ5-B
  SEQ6-B

Request a maximum of 5 strains. 4 strains should pass, even though this means groups are sampled unevenly.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --group-by group \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 0 \
  >  --output-strains filtered_strains.txt
  Sampling at 3 per group.
  2 strains were dropped during filtering
  	2 were dropped because of subsampling criteria
  4 strains passed all filters
  $ sort filtered_strains.txt
  SEQ1-A
  SEQ2-B
  SEQ5-B
  SEQ6-B

Request a maximum of 3 strains. 2 strains should pass, one from each group.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --group-by group \
  >  --subsample-max-sequences 3 \
  >  --subsample-seed 0 \
  >  --output-strains filtered_strains.txt
  Sampling at 1 per group.
  4 strains were dropped during filtering
  	4 were dropped because of subsampling criteria
  2 strains passed all filters
  $ sort filtered_strains.txt
  SEQ1-A
  SEQ6-B
