
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur ancestral --tree $TEST_DATA_DIR/in/tree.nwk --alignment $TEST_DATA_DIR/in/masked.vcf.gz --output-node-data $TMP/out/nt_muts.json --inference joint --output-vcf $TMP/out/nt_muts.vcf --vcf-reference $TEST_DATA_DIR/in/ref.fasta >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/nt_muts.json $TMP/out/nt_muts.json
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/nt_muts.vcf $TMP/out/nt_muts.vcf
  $ echo $?
  0