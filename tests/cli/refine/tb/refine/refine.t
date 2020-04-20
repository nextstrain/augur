
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur refine --tree $TEST_DATA_DIR/in/tree_raw.nwk --alignment $TEST_DATA_DIR/in/masked.vcf.gz --metadata $TEST_DATA_DIR/in/meta.tsv --output-tree $TMP/out/tree.nwk --output-node-data $TMP/out/branch_lengths.json --vcf-reference $TEST_DATA_DIR/in/ref.fasta --timetree --root min_dev --coalescent 0.0001 --clock-rate 1e-07 --clock-std-dev 3e-08 >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/tree.nwk $TMP/out/tree.nwk
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/branch_lengths.json $TMP/out/branch_lengths.json
  $ echo $?
  0