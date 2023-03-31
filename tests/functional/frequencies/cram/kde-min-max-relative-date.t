Setup

  $ source "$TESTDIR"/_setup.sh

Calculate KDE-based tip frequencies for a time period with relative dates.
Testing relative dates deterministically from the shell is tricky.
To keep these tests simple and avoid freezing the system date to specific values, this test checks the logical consistency of the requested relative dates and pivot interval.
With a minimum date of 12 months ago, a maximum date of 6 months ago, and a pivot interval of 3 months, we always expect 3 pivots in the final output corresponding to the end date (always included), the 3-month pivot before that end date, and the start date (also should always be included).
Since the test data are much older than the time period requested, all strains will always have frequencies of 0 for the 3 requested pivots.
As long as we always calculate 3 pivots with frequencies of 0 for all strains, we can ignore the actual pivot values calculated for the relative dates in the diff below.

  $ ${AUGUR} frequencies \
  >  --method kde \
  >  --tree "$TESTDIR/../data/tree.nwk" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --pivot-interval 3 \
  >  --min-date 12M \
  >  --max-date 6M \
  >  --output tip-frequencies.json > /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >  --exclude-paths "root['generated_by']['version']" "root['pivots']" -- \
  >  "$TESTDIR/../data/zika_tip-frequencies_with_relative_dates.json" \
  >  tip-frequencies.json
  {}
  $ rm -f tip-frequencies.json
