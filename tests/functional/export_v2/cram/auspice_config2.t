Setup

  $ source "$TESTDIR"/_setup.sh

Export with auspice config JSON with an extensions block

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --node-data "$TESTDIR/../data/div_node-data.json" "$TESTDIR/../data/location_node-data.json" \
  >   --auspice-config "$TESTDIR/../data/auspice_config2.json" \
  >   --output dataset2.json &>/dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py"  "$TESTDIR/../data/dataset2.json" dataset2.json \
  >   --exclude-paths "root['meta']['updated']"
  {}
