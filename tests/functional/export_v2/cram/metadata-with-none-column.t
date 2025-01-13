Setup

  $ source "$TESTDIR"/_setup.sh

Run export with metadata that contains "none" column and asked to use as coloring.
This is expected to output a warning that "none" is an invalid coloring column and skip it in colorings.

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --node-data "$TESTDIR/../data/div_node-data.json" "$TESTDIR/../data/location_node-data.json" \
  >   --metadata "$TESTDIR/../data/none_column_metadata.tsv" \
  >   --color-by-metadata "none" \
  >   --maintainers "Nextstrain Team" \
  >   --output dataset1.json >/dev/null
  WARNING: You asked for a color-by for trait 'none', but this is an invalid coloring key.
  It will be ignored during export, please rename field if you would like to use it as a coloring.
  \s{0} (re)

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-without-none-column.json" dataset1.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['maintainers']"
  {}

Run export with metadata that contains "none" column and asked to use as metadata.
This is expected to output a warning that "none" is an invalid node_attr and skip it in metadata.

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --node-data "$TESTDIR/../data/div_node-data.json" "$TESTDIR/../data/location_node-data.json" \
  >   --metadata "$TESTDIR/../data/none_column_metadata.tsv" \
  >   --metadata-column "none" \
  >   --maintainers "Nextstrain Team" \
  >   --output dataset2.json >/dev/null
  WARNING: You asked for a metadata field 'none', but this is an invalid field.
  It will be ignored during export, please rename field if you would like to include as a metadata field.
  \s{0} (re)

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-without-none-column.json" dataset2.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['maintainers']"
  {}
