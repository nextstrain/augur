Integration tests for augur frequencies.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="../../bin/augur"

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
With a minimum date of 1 year (12 months) ago, a maximum date of 6 months ago, and a pivot interval of 3 months, we expect only 2 pivots in the final output.

  $ ${AUGUR} frequencies \
  >  --method kde \
  >  --tree "frequencies/tree.nwk" \
  >  --metadata "frequencies/metadata.tsv" \
  >  --pivot-interval 3 \
  >  --min-date 1Y \
  >  --max-date 6M \
  >  --output "$TMP/tip-frequencies.json" > /dev/null

Since the test data are much older than the time period requested, all strains will have frequencies of 0 for the 2 requested pivots.
We can ignore the values of the calculated pivots which will vary based on when the test is run.

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" \
  >  --exclude-paths "root['generated_by']['version']" "root['pivots']" -- \
  >  "frequencies/zika_tip-frequencies_with_relative_dates.json" \
  >  "$TMP/tip-frequencies.json"
  {}
  $ rm -f "$TMP/tip-frequencies.json"

  $ popd > /dev/null
