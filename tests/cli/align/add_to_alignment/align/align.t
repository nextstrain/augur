
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur align --existing-alignment $TEST_DATA_DIR/in/aligned.fasta --sequences $TEST_DATA_DIR/in/new1.fasta $TEST_DATA_DIR/in/new2.fasta $TEST_DATA_DIR/in/new3_present_in_alignment.fasta $TEST_DATA_DIR/in/new4_present_in_seqs.fasta --reference-sequence $TEST_DATA_DIR/in/reference.gb --output $TMP/out/aligned.fasta --fill-gaps >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/aligned.fasta $TMP/out/aligned.fasta
  $ echo $?
  0