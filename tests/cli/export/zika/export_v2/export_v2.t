
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur export v2 --tree $TEST_DATA_DIR/in/tree.nwk --metadata $TEST_DATA_DIR/in/metadata.tsv --node-data $TEST_DATA_DIR/in/branch_lengths.json $TEST_DATA_DIR/in/traits.json $TEST_DATA_DIR/in/nt_muts.json $TEST_DATA_DIR/in/aa_muts.json --colors $TEST_DATA_DIR/config/colors.tsv --auspice-config $TEST_DATA_DIR/config/auspice_config_v2.json --output $TMP/out/v2_zika.json --title 'Real-time tracking of Zika virus evolution -- v2 JSON' --panels tree map entropy frequencies >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/v2_zika.json $TMP/out/v2_zika.json
  $ echo $?
  0