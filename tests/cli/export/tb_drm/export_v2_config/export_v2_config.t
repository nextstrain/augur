
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur export v2 --tree $TEST_DATA_DIR/in/tree.nwk --metadata $TEST_DATA_DIR/in/meta.tsv --node-data $TEST_DATA_DIR/in/branch_lengths.json $TEST_DATA_DIR/in/drms.json $TEST_DATA_DIR/in/traits.json $TEST_DATA_DIR/in/aa_muts.json $TEST_DATA_DIR/in/nt_muts.json --colors $TEST_DATA_DIR/in/color.tsv --lat-longs $TEST_DATA_DIR/in/lat_longs.tsv --output $TMP/out/v2_tbdrm_config.json --auspice-config $TEST_DATA_DIR/in/v2_config.json >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/v2_tbdrm_config.json $TMP/out/v2_tbdrm_config.json
  $ echo $?
  0