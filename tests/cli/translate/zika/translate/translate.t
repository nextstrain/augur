
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

Explicitly cd into TEST_DATA_DIR, so we can pass a relative path to --reference-sequence.  
The value for this param gets written into the output file.  So if we just passed the full path, e.g. $TEST_DATA_DIR/in/foobar.xyz,
this test would fail on other machines/environments.

  $ cd $TEST_DATA_DIR
  $ augur translate --tree $TEST_DATA_DIR/in/tree.nwk --ancestral-sequences $TEST_DATA_DIR/in/nt_muts.json --reference-sequence config/zika_outgroup.gb --output-node-data $TMP/out/aa_muts.json >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/aa_muts.json $TMP/out/aa_muts.json
  $ echo $?
  0
