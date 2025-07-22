Setup

  $ source "$TESTDIR"/_setup.sh

Use the same options with 3 different compression methods.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 0 \
  >  --output-metadata filtered_metadata.tsv.gz 2>/dev/null

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 0 \
  >  --output-metadata filtered_metadata.tsv.xz 2>/dev/null

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 0 \
  >  --output-metadata filtered_metadata.tsv.zst 2>/dev/null

# The uncompressed outputs are identical.

  $ diff <(gzcat filtered_metadata.tsv.gz) <(xzcat filtered_metadata.tsv.xz)

  $ diff <(gzcat filtered_metadata.tsv.gz) <(zstdcat filtered_metadata.tsv.zst)
