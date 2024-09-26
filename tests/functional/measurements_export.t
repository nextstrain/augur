Integration tests for augur measurements export.

  $ source "$TESTDIR"/_setup.sh
  $ pushd "$TESTDIR" > /dev/null

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
  [1]

Minimal measurements export with user provided strain, value, and subset of columns.

  $ ${AUGUR} measurements export \
  >   --collection measurements_export/collection_without_strain_value_columns.tsv \
  >   --strain-column strain_field \
  >   --value-column value_field \
  >   --grouping-column field_1 \
  >   --include-columns field_1 field_3 \
  >   --output-json "$TMP/minimal_measurements_subset.json" &>/dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" measurements_export/minimal_measurements_subset.json "$TMP/minimal_measurements_subset.json"
  {}

Try measurements export with grouping column missing from include columns list
This is expected to fail.

  $ ${AUGUR} measurements export \
  >   --collection measurements_export/collection_without_strain_value_columns.tsv \
  >   --strain-column strain_field \
  >   --value-column value_field \
  >   --grouping-column field_1 \
  >   --include-columns field_3 \
  >   --output-json "$TMP/minimal_measurements_subset.json" 1>/dev/null
  ERROR: Provided grouping column 'field_1' was not in the list of columns to include: ['field_3'].
  [1]

Try measurements export with invalid grouping columns.
This is expected to fail.

  $ ${AUGUR} measurements export \
  >   --collection measurements_export/collection.tsv \
  >   --grouping-column bad_field \
  >   --output-json "$TMP/minimal_measurements.json"
  ERROR: Provided grouping column 'bad_field' does not exist in collection TSV.
  [1]

Measurements export for a single collection using only command line configs.

  $ ${AUGUR} measurements export \
  >   --collection measurements_export/collection.tsv \
  >   --grouping-column field_1 field_2 \
  >   --key args-collection \
  >   --title collection-display-title \
  >   --x-axis-label label \
  >   --threshold 2.0 \
  >   --filters field_1 field_2 \
  >   --group-by field_1 \
  >   --measurements-display mean \
  >   --show-overall-mean \
  >   --show-threshold \
  >   --output-json "$TMP/single_collection_with_args_measurements.json" &>/dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" measurements_export/single_collection_with_args_measurements.json "$TMP/single_collection_with_args_measurements.json"
  {}

Measurements export for a single collection using a collection config.

  $ ${AUGUR} measurements export \
  >   --collection measurements_export/collection.tsv \
  >   --collection-config measurements_export/collection_config.json \
  >   --output-json "$TMP/single_collection_with_config_measurements.json" &>/dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" measurements_export/single_collection_with_config_measurements.json "$TMP/single_collection_with_config_measurements.json"
  {}

Measurements export for a single collection using a collection config and command-line overrides.

  $ ${AUGUR} measurements export \
  >   --collection measurements_export/collection.tsv \
  >   --collection-config measurements_export/collection_config.json \
  >   --grouping-column field_3 \
  >   --key override-collection \
  >   --title override-collection-display-title \
  >   --x-axis-label override-label \
  >   --threshold 10.0 \
  >   --filters field_3 \
  >   --group-by field_3 \
  >   --measurements-display raw \
  >   --hide-overall-mean \
  >   --hide-threshold \
  >   --output-json "$TMP/single_collection_with_overrides_measurements.json" &>/dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" measurements_export/single_collection_with_overrides_measurements.json "$TMP/single_collection_with_overrides_measurements.json"
  {}

Measurements export for a single collection using a collection config and command-line overrides with multiple thresholds.

  $ ${AUGUR} measurements export \
  >   --collection measurements_export/collection.tsv \
  >   --collection-config measurements_export/collection_config.json \
  >   --grouping-column field_3 \
  >   --key override-collection \
  >   --title override-collection-display-title \
  >   --x-axis-label override-label \
  >   --thresholds 2.0 0.0 \
  >   --filters field_3 \
  >   --group-by field_3 \
  >   --measurements-display raw \
  >   --hide-overall-mean \
  >   --hide-threshold \
  >   --output-json "$TMP/single_collection_with_multiple_thresholds.json" &>/dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" measurements_export/single_collection_with_multiple_thresholds.json "$TMP/single_collection_with_multiple_thresholds.json"
  {}
