
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur tree --alignment $TEST_DATA_DIR/in/aligned.fasta --output $TMP/out/tree_raw.nwk --method iqtree >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/tree_raw.nwk $TMP/out/tree_raw.nwk
  $ echo $?
  0