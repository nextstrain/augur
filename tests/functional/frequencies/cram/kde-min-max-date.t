Setup

  $ source "$TESTDIR"/_setup.sh

Calculate KDE-based tip frequencies for a time period with fixed dates.
Pivots get calculated from the requested date range.

  $ ${AUGUR} frequencies \
  >  --method kde \
  >  --tree "$TESTDIR/../data/tree.nwk" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --pivot-interval 3 \
  >  --min-date 2015-01-01 \
  >  --max-date 2016-01-01 \
  >  --output tip-frequencies.json > /dev/null

  $ diff -u --ignore-matching-lines version "$TESTDIR/../data/zika_tip-frequencies_with_fixed_dates.json" tip-frequencies.json
  $ rm -f tip-frequencies.json
