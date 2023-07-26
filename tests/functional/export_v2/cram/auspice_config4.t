Setup

  $ source "$TESTDIR"/_setup.sh

Run export with metadata and external colors TSV that contains zero values.

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --node-data "$TESTDIR/../data/div_node-data.json" "$TESTDIR/../data/location_node-data.json" \
  >   --auspice-config "$TESTDIR/../data/auspice_config4.json" \
  >   --metadata "$TESTDIR/../data/zero_value_metadata.tsv" \
  >   --colors "$TESTDIR/../data/zero_value_colors.tsv" \
  >   --output dataset.json &> /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset3.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['maintainers']"
  {}
