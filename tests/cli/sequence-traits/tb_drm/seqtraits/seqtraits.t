
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur sequence-traits --ancestral-sequences $TEST_DATA_DIR/in/nt_muts.vcf --vcf-reference $TEST_DATA_DIR/in/ref.fasta --translations $TEST_DATA_DIR/in/translations.vcf --vcf-translate-reference $TEST_DATA_DIR/in/translations_reference.fasta --features $TEST_DATA_DIR/in/DRMs-AAnuc.tsv --output-node-data $TMP/out/drms.json --count traits --label Drug_Resistance >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/drms.json $TMP/out/drms.json
  $ echo $?
  0