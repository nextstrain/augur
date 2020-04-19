
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur export v1 --tree $TEST_DATA_DIR/in/tree.nwk --metadata $TEST_DATA_DIR/in/meta.tsv --node-data $TEST_DATA_DIR/in/branch_lengths.json $TEST_DATA_DIR/in/drms.json $TEST_DATA_DIR/in/traits.json $TEST_DATA_DIR/in/aa_muts.json $TEST_DATA_DIR/in/nt_muts.json --auspice-config $TEST_DATA_DIR/in/config.json --colors $TEST_DATA_DIR/in/color.tsv --lat-longs $TEST_DATA_DIR/in/lat_longs.tsv --output-tree $TMP/out/v1_tbdrm_tree.json --output-meta $TMP/out/v1_tbdrm_meta.json >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/v1_tbdrm_tree.json $TMP/out/v1_tbdrm_tree.json
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/v1_tbdrm_meta.json $TMP/out/v1_tbdrm_meta.json
  $ echo $?
  0