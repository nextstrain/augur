
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur ancestral --tree $TEST_DATA_DIR/in/tree.nwk --alignment $TEST_DATA_DIR/in/aligned.fasta --output-node-data $TMP/out/nt_muts.json --inference joint >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/nt_muts.json $TMP/out/nt_muts.json
  $ echo $?
  0