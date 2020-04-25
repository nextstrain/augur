
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur export v2 --tree $TEST_DATA_DIR/in/tree_temporal.nwk --node-data $TEST_DATA_DIR/in/time_only_branch_lengths.json --metadata $TEST_DATA_DIR/in/metadata_boolean.tsv --color-by-metadata top_half_True_False top_half_1_0 top_half_yes_no all_missing --colors $TEST_DATA_DIR/../zika/config/colors.tsv --title "Boolean metadata traits" --output $TMP/out/v2_boolean-metadata.json >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/v2_boolean-metadata.json $TMP/out/v2_boolean-metadata.json
  $ echo $?
  0