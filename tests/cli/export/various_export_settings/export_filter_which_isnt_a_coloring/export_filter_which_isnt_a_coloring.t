
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur export v2 --tree $TEST_DATA_DIR/in/tree_temporal.nwk --node-data $TEST_DATA_DIR/in/temporal_branch_lengths.json --metadata $TEST_DATA_DIR/in/metadata.tsv --auspice-config $TEST_DATA_DIR/config/filters.json --output $TMP/out/v2_filters-not-colors.json >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/v2_filters-not-colors.json $TMP/out/v2_filters-not-colors.json
  $ echo $?
  0