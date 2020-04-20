
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

Explicitly cd into TEST_DATA_DIR, so we can pass a relative path to --tree AND --alignment.  
The value for this param gets written into the output file.  So if we just passed the full path, e.g. $TEST_DATA_DIR/in/tree_raw.nwk,
this test would fail on other machines/environments.

  $ cd $TEST_DATA_DIR
  $ augur refine --tree in/tree_raw.nwk --alignment in/aligned.fasta --metadata $TEST_DATA_DIR/in/metadata.tsv --output-tree $TMP/out/tree_temporal.nwk --output-node-data $TMP/out/temporal_branch_lengths.json --timetree --coalescent opt --date-confidence --date-inference marginal --clock-filter-iqd 4 >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/tree_temporal.nwk $TMP/out/tree_temporal.nwk
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/temporal_branch_lengths.json $TMP/out/temporal_branch_lengths.json
  $ echo $?
  0
