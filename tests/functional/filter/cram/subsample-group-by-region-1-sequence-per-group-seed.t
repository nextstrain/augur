Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Pandas engine
-------------

Filter with subsampling, requesting 1 sequence per group (for a group with 4 distinct values).

  $ ${AUGUR} filter \
  >  --metadata filter/data/metadata.tsv \
  >  --group-by region \
  >  --sequences-per-group 1 \
  >  --subsample-seed 314159 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  $ wc -l "$TMP/filtered_strains.txt"
  \s*4 .* (re)

By setting the subsample seed above, we should guarantee that we get the same "random" strains as another run with the same command.

  $ ${AUGUR} filter \
  >  --metadata filter/data/metadata.tsv \
  >  --group-by region \
  >  --sequences-per-group 1 \
  >  --subsample-seed 314159 \
  >  --output-strains "$TMP/filtered_strains_repeated.txt" > /dev/null

  $ diff -u <(sort "$TMP/filtered_strains.txt") <(sort "$TMP/filtered_strains_repeated.txt")
  $ rm -f "$TMP/filtered_strains.txt" "$TMP/filtered_strains_repeated.txt"

SQLite engine
-------------

Filter with subsampling, requesting 1 sequence per group (for a group with 4 distinct values).

  $ ${AUGUR} filter --engine sqlite \
  >  --metadata filter/data/metadata.tsv \
  >  --group-by region \
  >  --sequences-per-group 1 \
  >  --subsample-seed 314159 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  $ wc -l "$TMP/filtered_strains.txt"
  \s*4 .* (re)

By setting the subsample seed above, we should guarantee that we get the same "random" strains as another run with the same command.

  $ ${AUGUR} filter --engine sqlite \
  >  --metadata filter/data/metadata.tsv \
  >  --group-by region \
  >  --sequences-per-group 1 \
  >  --subsample-seed 314159 \
  >  --output-strains "$TMP/filtered_strains_repeated.txt" > /dev/null

  $ diff -u <(sort "$TMP/filtered_strains.txt") <(sort "$TMP/filtered_strains_repeated.txt")
  $ rm -f "$TMP/filtered_strains.txt" "$TMP/filtered_strains_repeated.txt"
