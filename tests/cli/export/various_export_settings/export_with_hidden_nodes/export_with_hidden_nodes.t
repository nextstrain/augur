
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur export v2 --tree $TEST_DATA_DIR/in/tree_temporal.nwk --node-data $TEST_DATA_DIR/in/temporal_branch_lengths.json $TEST_DATA_DIR/in/hidden_nodes.json --metadata $TEST_DATA_DIR/in/metadata.tsv --colors $TEST_DATA_DIR/../zika/config/colors.tsv --auspice-config $TEST_DATA_DIR/../zika/config/auspice_config_v2.json --title "Tree with hidden nodes (both div + temporal)" --output $TMP/out/v2_hidden-nodes.json >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/v2_hidden-nodes.json $TMP/out/v2_hidden-nodes.json
  $ echo $?
  0