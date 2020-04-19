
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur export v2 --tree $TEST_DATA_DIR/in/tree.nwk --metadata $TEST_DATA_DIR/in/metadata.tsv --node-data $TEST_DATA_DIR/in/branch_lengths.json $TEST_DATA_DIR/in/nt_muts.json $TEST_DATA_DIR/in/aa_muts.json --auspice-config $TEST_DATA_DIR/in/auspice_config.json --output $TMP/out/progressive_align.json >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/progressive_align.json $TMP/out/progressive_align.json
  $ echo $?
  0