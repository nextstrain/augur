
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur filter --sequences $TEST_DATA_DIR/in/drm.vcf.gz --metadata $TEST_DATA_DIR/in/meta.tsv --output $TMP/out/filtered.vcf.gz --exclude $TEST_DATA_DIR/in/dropped_strains.txt >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/filtered.vcf.gz $TMP/out/filtered.vcf.gz
  $ echo $?
  0