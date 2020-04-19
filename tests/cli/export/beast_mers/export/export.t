
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur export v2 --tree $TEST_DATA_DIR/in/tree.new --node-data $TEST_DATA_DIR/in/beast_data.json --auspice-config $TEST_DATA_DIR/config/auspice_config.json --output $TMP/out/v2_mers-cov.json >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/v2_mers-cov.json $TMP/out/v2_mers-cov.json
  $ echo $?
  0