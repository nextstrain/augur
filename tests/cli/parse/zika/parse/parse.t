
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur parse --sequences $TEST_DATA_DIR/in/zika.fasta --output-sequences $TMP/out/sequences.fasta --output-metadata $TMP/out/metadata.tsv --fields strain virus accession date region country division city db segment authors url title journal paper_url --prettify-fields region country division city >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/sequences.fasta $TMP/out/sequences.fasta
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/metadata.tsv $TMP/out/metadata.tsv
  $ echo $?
  0