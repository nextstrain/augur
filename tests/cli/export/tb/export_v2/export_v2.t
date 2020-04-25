
  $ TEST_DATA_DIR="$TESTDIR"
  $ mkdir -p "$TMP/out"

  $ augur export v2 --tree $TEST_DATA_DIR/in/tree.nwk --metadata $TEST_DATA_DIR/in/meta.tsv --node-data $TEST_DATA_DIR/in/branch_lengths.json $TEST_DATA_DIR/in/traits.json $TEST_DATA_DIR/in/aa_muts.json $TEST_DATA_DIR/in/nt_muts.json $TEST_DATA_DIR/in/clades.json --colors $TEST_DATA_DIR/in/color.tsv --lat-longs $TEST_DATA_DIR/in/lat_longs.tsv --output $TMP/out/tb_main.json --title 'TB outbreak in Nunavik, Canada' --maintainers 'Emma Hodcroft <https://neherlab.org/emma-hodcroft.html>' 'John Brown <http://www.google.com>' --geo-resolutions location >/dev/null
  $ echo $?
  0
  $ diff -u $TEST_DATA_DIR/expected/tb_main.json $TMP/out/tb_main.json
  $ echo $?
  0