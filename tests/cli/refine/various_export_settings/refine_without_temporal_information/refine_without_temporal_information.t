
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

Explicitly cd into TEST_DATA_DIR, so we can pass a relative path to --tree AND --alignment.  
The value for this param gets written into the output file.  So if we just passed the full path, e.g. $TEST_DATA_DIR/in/tree_raw.nwk,
this test would fail on other machines/environments.

  $ cd $TEST_DATA_DIR
  $ augur refine --tree in/tree_raw.nwk --alignment in/aligned.fasta --metadata $TEST_DATA_DIR/in/metadata.tsv --output-tree $TMP/out/tree_div_only.nwk --output-node-data $TMP/out/div_only_branch_lengths.json --keep-root >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/tree_div_only.nwk $TMP/out/tree_div_only.nwk
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/div_only_branch_lengths.json $TMP/out/div_only_branch_lengths.json
  $ echo $?
  0
