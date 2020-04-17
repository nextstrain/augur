
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur traits --tree $TEST_DATA_DIR/in/tree.nwk --metadata $TEST_DATA_DIR/in/meta.tsv --columns country region --output-node-data $TMP/out/traits.json >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/traits.json $TMP/out/traits.json
  $ echo $?
  0