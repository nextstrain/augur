Integration tests for augur measurements export.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="../../bin/augur"

Minimal measurements export with existing strain and value columns.

  $ ${AUGUR} measurements export \
  >   --collection measurements_export/collection.tsv \
  >   --grouping-column field_1 \
  >   --output-json "$TMP/minimal_measurements.json" &>/dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" measurements_export/minimal_measurements.json "$TMP/minimal_measurements.json"
  {}

Minimal measurements export with user provided strain and value columns.

  $ ${AUGUR} measurements export \
  >   --collection measurements_export/collection_without_strain_value_columns.tsv \
  >   --strain-column strain_field \
  >   --value-column value_field \
  >   --grouping-column field_1 \
  >   --output-json "$TMP/minimal_measurements.json" &>/dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" measurements_export/minimal_measurements.json "$TMP/minimal_measurements.json" \
  >   --exclude-paths "root['collections'][0]['key']"
  {}

Try measurements export with user provided strain and value columns that would overwrite existing columns.
This is expected to fail.

  $ ${AUGUR} measurements export \
  >   --collection measurements_export/collection.tsv \
  >   --strain-column field_1 \
  >   --value-column field_2 \
  >   --grouping-column field_1 \
  >   --output-json "$TMP/minimal_measurements.json"
  ERROR: Cannot use provided 'field_1' column as the strain column because a 'strain' column already exists in collection TSV.
  ERROR: Cannot use provided 'field_2' column as the value column because a 'value' column already exists in collection TSV.
  ERROR: Loading of collection TSV was unsuccessful. See detailed errors above.
  [1]

Try measurements export with user provided strain and value columns that are the same column.
This is expected to fail.

  $ ${AUGUR} measurements export \
  >   --collection measurements_export/collection_without_strain_value_columns.tsv \
  >   --strain-column field_1 \
  >   --value-column field_1 \
  >   --grouping-column field_1 \
  >   --output-json "$TMP/minimal_measurements.json"
  ERROR: The strain column and value column cannot be the same column.
  ERROR: Loading of collection TSV was unsuccessful. See detailed errors above.
  [1]
