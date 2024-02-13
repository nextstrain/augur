Setup

  $ source "$TESTDIR"/_setup.sh

Calculate diffusion-based tip frequencies from a refined tree.
This currently doesn't work due to a bug.

  $ ${AUGUR} frequencies \
  >  --method diffusion \
  >  --tree "$TESTDIR/../data/tree.nwk" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --pivot-interval 3 \
  >  --output tip-frequencies.json > /dev/null 2>&1
  [2]
