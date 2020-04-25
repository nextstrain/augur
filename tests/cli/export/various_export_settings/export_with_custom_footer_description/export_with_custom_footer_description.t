
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur export v2 --tree $TEST_DATA_DIR/in/tree_temporal.nwk --node-data $TEST_DATA_DIR/in/temporal_branch_lengths.json --metadata $TEST_DATA_DIR/in/metadata.tsv --description $TEST_DATA_DIR/config/footer-description.md --title "Custom footer" --output $TMP/out/v2_custom-footer.json >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/v2_custom-footer.json $TMP/out/v2_custom-footer.json
  $ echo $?
  0