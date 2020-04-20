
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur filter --sequences $TEST_DATA_DIR/in/sequences.fasta --metadata $TEST_DATA_DIR/in/metadata.tsv --exclude $TEST_DATA_DIR/config/dropped_strains.txt --output $TMP/out/filtered.fasta --group-by country year month --sequences-per-group 1 --min-date 2012 >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/filtered.fasta $TMP/out/filtered.fasta
  $ echo $?
  0