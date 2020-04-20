
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur refine --tree $TEST_DATA_DIR/in/tree_raw.nwk --alignment $TEST_DATA_DIR/in/aligned.fasta --metadata $TEST_DATA_DIR/in/metadata.tsv --output-tree $TMP/out/tree_div_only.nwk --output-node-data $TMP/out/div_only_branch_lengths.json --keep-root >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/tree_div_only.nwk $TMP/out/tree_div_only.nwk
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/div_only_branch_lengths.json $TMP/out/div_only_branch_lengths.json
  $ echo $?
  0