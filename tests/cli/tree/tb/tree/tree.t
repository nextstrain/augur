
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur tree --exclude-sites $TEST_DATA_DIR/in/drm_sites.txt --alignment $TEST_DATA_DIR/in/masked.vcf.gz --vcf-reference $TEST_DATA_DIR/in/ref.fasta --output $TMP/out/tree_raw.nwk --method fasttree >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/tree_raw.nwk $TMP/out/tree_raw.nwk
  $ echo $?
  0