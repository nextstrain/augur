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