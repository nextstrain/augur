
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

Explicitly cd into TEST_DATA_DIR, so we can pass a relative path to --mcc.  
The value for this param gets written into the output file.  So if we instead passed the full path as usual, e.g. $TEST_DATA_DIR/in/foobar.tree,
this test would fail on other machines/environments.

  $ cd $TEST_DATA_DIR
  $ augur import beast --mcc in/MERS_CoV_274_mcc.tree --output-tree $TMP/out/tree.new --output-node-data $TMP/out/beast_data.json >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/beast_data.json $TMP/out/beast_data.json
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/tree.new $TMP/out/tree.new
  $ echo $?
  0
