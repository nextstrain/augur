
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur refine --tree $TEST_DATA_DIR/in/tree_raw.nwk --alignment $TEST_DATA_DIR/in/aligned.fasta --metadata $TEST_DATA_DIR/in/metadata.tsv --output-tree $TMP/out/tree_temporal.nwk --output-node-data $TMP/out/temporal_branch_lengths.json --timetree --coalescent opt --date-confidence --date-inference marginal --clock-filter-iqd 4 >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/tree_temporal.nwk $TMP/out/tree_temporal.nwk
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/temporal_branch_lengths.json $TMP/out/temporal_branch_lengths.json
  $ echo $?
  0