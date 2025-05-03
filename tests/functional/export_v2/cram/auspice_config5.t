Setup

  $ source "$TESTDIR"/_setup.sh

Run export with metadata and auspice config that include deprecated field names.

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --node-data "$TESTDIR/../data/div_node-data.json" \
  >   --auspice-config "$TESTDIR/../data/auspice_config5.json" \
  >   --metadata "$TESTDIR/../data/deprecated_metadata.tsv" \
  >   --output dataset.json &> /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset5.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['maintainers']"
  {}
