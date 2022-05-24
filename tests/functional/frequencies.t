Integration tests for augur frequencies.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="${AUGUR:-../../bin/augur}"

Calculate KDE-based tip frequencies from a refined tree.
Timepoints used to estimate frequencies (i.e., "pivots") get calculated from the range of dates in the given metadata.

  $ ${AUGUR} frequencies \
  >  --method kde \
  >  --tree "frequencies/tree.nwk" \
  >  --metadata "frequencies/metadata.tsv" \
  >  --pivot-interval 3 \
  >  --output "$TMP/tip-frequencies.json" > /dev/null

  $ diff -u --ignore-matching-lines version "frequencies/zika_tip-frequencies.json" "$TMP/tip-frequencies.json"
  $ rm -f "$TMP/tip-frequencies.json"

Calculate KDE-based tip frequencies for a time period with fixed dates.
Pivots get calculated from the requested date range.

  $ ${AUGUR} frequencies \
  >  --method kde \
  >  --tree "frequencies/tree.nwk" \
  >  --metadata "frequencies/metadata.tsv" \
  >  --pivot-interval 3 \
  >  --min-date 2015-01-01 \
  >  --max-date 2016-01-01 \
  >  --output "$TMP/tip-frequencies.json" > /dev/null

  $ diff -u --ignore-matching-lines version "frequencies/zika_tip-frequencies_with_fixed_dates.json" "$TMP/tip-frequencies.json"
  $ rm -f "$TMP/tip-frequencies.json"

Calculate KDE-based tip frequencies for a time period with relative dates.
Testing relative dates deterministically from the shell is tricky.
To keep these tests simple and avoid freezing the system date to specific values, this test checks the logical consistency of the requested relative dates and pivot interval.
With a minimum date of 1 year (12 months) ago, a maximum date of 6 months ago, and a pivot interval of 3 months, we expect only 2 pivots in the final output.
Since the test data are much older than the time period requested, all strains will always have frequencies of 0 for the 2 requested pivots.
As long as we always calculate 2 pivots with frequencies of 0 for all strains, we can ignore the actual pivot values calculated for the relative dates in the diff below.

  $ ${AUGUR} frequencies \
  >  --method kde \
  >  --tree "frequencies/tree.nwk" \
  >  --metadata "frequencies/metadata.tsv" \
  >  --pivot-interval 3 \
  >  --min-date 1Y \
  >  --max-date 6M \
  >  --output "$TMP/tip-frequencies.json" > /dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" \
  >  --exclude-paths "root['generated_by']['version']" "root['pivots']" -- \
  >  "frequencies/zika_tip-frequencies_with_relative_dates.json" \
  >  "$TMP/tip-frequencies.json"
  {}
  $ rm -f "$TMP/tip-frequencies.json"

  $ popd > /dev/null
