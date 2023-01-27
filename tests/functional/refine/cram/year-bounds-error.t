Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata with a strain that has partial ambiguity on the century-level
(20XX) so that --year-bounds is applied.

  $ cat >metadata.tsv <<~~
  > strain	date
  > PAN/CDC_259359_V1_V3/2015	20XX-XX-XX
  > ~~


Check that invalid --year-bounds provides useful error messages.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --metadata metadata.tsv \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --timetree \
  >  --year-bounds 1950 1960 1970 \
  >  --divergence-units mutations > /dev/null
  ERROR: Invalid value for --year-bounds: The year bounds [1950, 1960, 1970] must have only one (lower) or two (lower, upper) bounds.
  [2]

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --metadata metadata.tsv \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --timetree \
  >  --year-bounds 0 1000 \
  >  --divergence-units mutations > /dev/null
  ERROR: Invalid value for --year-bounds: 0 is not a valid year.
  [2]

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --metadata metadata.tsv \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --timetree \
  >  --year-bounds -2000 -3000 \
  >  --divergence-units mutations > /dev/null
  ERROR: Invalid value for --year-bounds: -3000 is not a valid year.
  [2]
