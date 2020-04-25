
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur export v2 --tree $TEST_DATA_DIR/in/tree_temporal.nwk --node-data $TEST_DATA_DIR/in/temporal_branch_lengths.json --metadata $TEST_DATA_DIR/in/metadata.tsv --colors $TEST_DATA_DIR/../zika/config/colors.tsv --auspice-config $TEST_DATA_DIR/config/default-layout.json --output $TMP/out/v2_default-layout.json >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/v2_default-layout.json $TMP/out/v2_default-layout.json
  $ echo $?
  0