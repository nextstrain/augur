Setup

  $ source "$TESTDIR"/_setup.sh

Calculate KDE-based tip frequencies from a refined tree.
Timepoints used to estimate frequencies (i.e., "pivots") get calculated from the range of dates in the given metadata.

  $ ${AUGUR} frequencies \
  >  --method kde \
  >  --tree "$TESTDIR/../data/tree.nwk" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --pivot-interval 3 \
  >  --output tip-frequencies.json > /dev/null

  $ diff -u --ignore-matching-lines version "$TESTDIR/../data/zika_tip-frequencies.json" tip-frequencies.json
  $ rm -f tip-frequencies.json
