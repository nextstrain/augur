
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur align --sequences $TEST_DATA_DIR/in/filtered.fasta --reference-sequence $TEST_DATA_DIR/config/zika_outgroup.gb --output $TMP/out/aligned.fasta --fill-gaps >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/aligned.fasta $TMP/out/aligned.fasta
  $ echo $?
  0