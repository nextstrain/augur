
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur mask --sequences $TEST_DATA_DIR/in/filtered.vcf.gz --output $TMP/out/masked.vcf.gz --mask $TEST_DATA_DIR/in/Locus_to_exclude_Mtb.bed >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/masked.vcf.gz $TMP/out/masked.vcf.gz
  $ echo $?
  0