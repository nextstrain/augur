
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur frequencies --method kde --tree $TEST_DATA_DIR/in/tree.nwk --metadata $TEST_DATA_DIR/in/metadata.tsv --pivot-interval 3 --output $TMP/out/zika_tip-frequencies.json >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/zika_tip-frequencies.json $TMP/out/zika_tip-frequencies.json
  $ echo $?
  0