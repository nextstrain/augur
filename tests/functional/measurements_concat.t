Integration tests for augur measurements export.

  $ source "$TESTDIR"/_setup.sh
  $ pushd "$TESTDIR" > /dev/null

Measurements concat for two measurements JSONs, each with a single collection.

  $ ${AUGUR} measurements concat \
  >   --jsons measurements_concat/single_collection_measurements_1.json measurements_concat/single_collection_measurements_2.json \
  >   --default-collection collection_1 \
  >   --output-json "$TMP/two_collections_measurements.json" &>/dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" measurements_concat/two_collections_measurements.json "$TMP/two_collections_measurements.json"
  {}

Measurements concat for two measurements JSONs, where one has multiple collections.

  $ ${AUGUR} measurements concat \
  >   --jsons measurements_concat/two_collections_measurements.json measurements_concat/single_collection_measurements_3.json \
  >   --default-collection collection_1 \
  >   --output-json "$TMP/multiple_collections_measurements.json" &>/dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" measurements_concat/multiple_collections_measurements.json "$TMP/multiple_collections_measurements.json"
  {}

Measurements concat for measurements JSONs that have collections that share the same key.
This is expected to fail.

  $ ${AUGUR} measurements concat \
  >   --jsons measurements_concat/single_collection_measurements_1.json measurements_concat/single_collection_measurements_1.json \
  >   --default-collection collection_1 \
  >   --output-json "$TMP/multiple_collections_measurements.json" 1>/dev/null
  ERROR: Collections at indexes [0, 1] share the same collection key 'collection_1'.
  ERROR: Validation of output JSON failed. See detailed errors above.
  [1]

Measurements concat with an invalid default collection.
This is expected to fail.

  $ ${AUGUR} measurements concat \
  >   --jsons measurements_concat/single_collection_measurements_1.json measurements_concat/single_collection_measurements_2.json \
  >   --default-collection collection_3 \
  >   --output-json "$TMP/multiple_collections_measurements.json" 1>/dev/null
  ERROR: The default collection key 'collection_3' does not match any of the collections' keys.
  ERROR: Validation of output JSON failed. See detailed errors above.
  [1]
