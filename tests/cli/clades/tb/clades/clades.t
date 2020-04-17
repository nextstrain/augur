
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur clades --tree $TEST_DATA_DIR/in/tree.nwk --mutations $TEST_DATA_DIR/in/nt_muts.json $TEST_DATA_DIR/in/aa_muts.json --output-node-data $TMP/out/clades.json --clades $TEST_DATA_DIR/in/clades.tsv >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/clades.json $TMP/out/clades.json
  $ echo $?
  0