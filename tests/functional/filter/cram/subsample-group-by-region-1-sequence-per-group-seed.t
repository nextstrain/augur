Setup

  $ source "$TESTDIR"/_setup.sh

Filter with subsampling, requesting 1 sequence per group (for a group with 4 distinct values).

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --group-by region \
  >  --sequences-per-group 1 \
  >  --subsample-seed 314159 \
  >  --output-strains filtered_strains.txt 2>/dev/null
  $ wc -l filtered_strains.txt
  \s*4 .* (re)

By setting the subsample seed above, we should guarantee that we get the same "random" strains as another run with the same command.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --group-by region \
  >  --sequences-per-group 1 \
  >  --subsample-seed 314159 \
  >  --output-strains filtered_strains_repeated.txt 2>/dev/null

  $ diff -u <(sort filtered_strains.txt) <(sort filtered_strains_repeated.txt)
