Setup

  $ source "$TESTDIR"/_setup.sh

Run export with metadata using the default id column of "strain".

  $ ${AUGUR} export v2 \
  >  --tree "$TESTDIR/../data/tree.nwk" \
  >  --metadata "$TESTDIR/../data/dataset1_metadata_with_strain.tsv" \
  >  --node-data "$TESTDIR/../data/div_node-data.json" "$TESTDIR/../data/location_node-data.json" \
  >  --auspice-config "$TESTDIR/../data/auspice_config1.json" \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset.json > /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset1.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['maintainers']"
  {}

Run export with metadata that uses a different id column other than "strain".
In this case, the column is "name" (one of the default columns expected by Augur's `io.read_metadata` function).

  $ ${AUGUR} export v2 \
  >  --tree "$TESTDIR/../data/tree.nwk" \
  >  --metadata "$TESTDIR/../data/dataset1_metadata_with_name.tsv" \
  >  --node-data "$TESTDIR/../data/div_node-data.json" "$TESTDIR/../data/location_node-data.json" \
  >  --auspice-config "$TESTDIR/../data/auspice_config1.json" \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset.json > /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset1.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['maintainers']"
  {}

Run export with metadata that uses an invalid id column.
This should fail with a helpful error message.

  $ ${AUGUR} export v2 \
  >  --tree "$TESTDIR/../data/tree.nwk" \
  >  --metadata "$TESTDIR/../data/dataset1_metadata_without_valid_id.tsv" \
  >  --node-data "$TESTDIR/../data/div_node-data.json" "$TESTDIR/../data/location_node-data.json" \
  >  --auspice-config "$TESTDIR/../data/auspice_config1.json" \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset.json > /dev/null
  ERROR: None of the possible id columns ('strain', 'name') were found in the metadata's columns ('invalid_id', 'div', 'mutation_length')
  [1]
