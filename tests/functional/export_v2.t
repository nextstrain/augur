Integration tests for augur export v2.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="../../bin/augur"

Minimal export

  $ ${AUGUR} export v2 \
  >   --tree export_v2/tree.nwk \
  >   --node-data export_v2/div_node-data.json \
  >   --output "$TMP/minimal.json" &>/dev/null
  [2]

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py"  export_v2/minimal.json "$TMP/minimal.json" \
  >   --exclude-paths "root['meta']['updated']"
  {}

Future test:
Run augur export _without_ any node-data JSONs when this can read divergence values from the newick file
and compare this to the tree using `div_node-data.json` - they should be identical.
See https://github.com/nextstrain/augur/issues/273 for more


Export with auspice config JSON which defines scale & legend settings
  $ ${AUGUR} export v2 \
  >   --tree export_v2/tree.nwk \
  >   --node-data export_v2/div_node-data.json export_v2/location_node-data.json \
  >   --auspice-config export_v2/auspice_config1.json \
  >   --output "$TMP/dataset1.json" &>/dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py"  export_v2/dataset1.json "$TMP/dataset1.json" \
  >   --exclude-paths "root['meta']['updated']"
  {}


Export with auspice config JSON with an extensions block
  $ ${AUGUR} export v2 \
  >   --tree export_v2/tree.nwk \
  >   --node-data export_v2/div_node-data.json export_v2/location_node-data.json \
  >   --auspice-config export_v2/auspice_config2.json \
  >   --output "$TMP/dataset2.json" &>/dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py"  export_v2/dataset2.json "$TMP/dataset2.json" \
  >   --exclude-paths "root['meta']['updated']"
  {}

# auspice_config3.json is the same as auspice_config2.json but with an extra key which the schema does not allow.
# Running without --skip-validation should result in an error
# Message printed: "Validation of export_v2/auspice_config3.json failed."
  $ ${AUGUR} export v2 \
  >   --tree export_v2/tree.nwk \
  >   --node-data export_v2/div_node-data.json export_v2/location_node-data.json \
  >   --auspice-config export_v2/auspice_config3.json \
  >   --output "$TMP/dataset2.json" &>/dev/null
  [2]

# Skipping validation gives us the same results as `auspice_config2.json`
  $ ${AUGUR} export v2 \
  >   --tree export_v2/tree.nwk \
  >   --node-data export_v2/div_node-data.json export_v2/location_node-data.json \
  >   --auspice-config export_v2/auspice_config3.json \
  >   --output "$TMP/dataset3.json" \
  >   --skip-validation &>/dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py"  export_v2/dataset2.json "$TMP/dataset3.json" \
  >   --exclude-paths "root['meta']['updated']"
  {}

Run export with metadata using the default id column of "strain".

  $ ${AUGUR} export v2 \
  >  --tree export_v2/tree.nwk \
  >  --metadata export_v2/dataset1_metadata_with_strain.tsv \
  >  --node-data export_v2/div_node-data.json export_v2/location_node-data.json \
  >  --auspice-config export_v2/auspice_config1.json \
  >  --maintainers "Nextstrain Team" \
  >  --output "$TMP/dataset1.json" > /dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" export_v2/dataset1.json "$TMP/dataset1.json" \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['maintainers']"
  {}
  $ rm -f "$TMP/dataset1.json"

Run export with metadata that uses a different id column other than "strain".
In this case, the column is "name" (one of the default columns expected by Augur's `io.read_metadata` function).

  $ ${AUGUR} export v2 \
  >  --tree export_v2/tree.nwk \
  >  --metadata export_v2/dataset1_metadata_with_name.tsv \
  >  --node-data export_v2/div_node-data.json export_v2/location_node-data.json \
  >  --auspice-config export_v2/auspice_config1.json \
  >  --maintainers "Nextstrain Team" \
  >  --output "$TMP/dataset1.json" > /dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" export_v2/dataset1.json "$TMP/dataset1.json" \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['maintainers']"
  {}
  $ rm -f "$TMP/dataset1.json"

Run export with metadata that uses an invalid id column.
This should fail with a helpful error message.

  $ ${AUGUR} export v2 \
  >  --tree export_v2/tree.nwk \
  >  --metadata export_v2/dataset1_metadata_without_valid_id.tsv \
  >  --node-data export_v2/div_node-data.json export_v2/location_node-data.json \
  >  --auspice-config export_v2/auspice_config1.json \
  >  --maintainers "Nextstrain Team" \
  >  --output "$TMP/dataset1.json" > /dev/null
  ERROR: None of the possible id columns (('strain', 'name')) were found in the metadata's columns ('invalid_id', 'div', 'mutation_length')
  [1]

  $ popd > /dev/null
