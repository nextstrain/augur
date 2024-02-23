Setup

  $ source "$TESTDIR"/_setup.sh

Calculate diffusion-based tip frequencies from a refined tree with `--regions`.
This currently does not work due to a bug.

  $ ${AUGUR} frequencies \
  >  --method diffusion \
  >  --tree "$TESTDIR/../data/tree.nwk" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --regions "global" "North America" "South America" \
  >  --pivot-interval 3 \
  >  --output tip-frequencies.json > /dev/null 2>&1
  [2]
