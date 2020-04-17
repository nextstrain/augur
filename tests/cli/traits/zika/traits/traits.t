
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur traits --tree $TEST_DATA_DIR/in/tree.nwk --weights config/trait_weights.csv --metadata $TEST_DATA_DIR/in/metadata.tsv --output-node-data $TMP/out/traits.json --columns country region --sampling-bias-correction 3 --confidence >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/traits.json $TMP/out/traits.json
  $ echo $?
  0