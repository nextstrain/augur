export AUGUR="${AUGUR:-$TESTDIR/../../../../bin/augur}"
set -o pipefail

# IQ-Tree writes to the input file directory, so we need to copy the data
# to allow running tests in parallel.
# This also avoids the data directory being polluted with output files.
cp -r "$TESTDIR/../data" "$TMP"
