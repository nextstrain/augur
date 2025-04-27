Setup

  $ source "$TESTDIR"/_setup.sh

Export with auspice config JSON which defines scale & legend settings

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --node-data "$TESTDIR/../data/div_node-data.json" "$TESTDIR/../data/location_node-data.json" \
  >   --auspice-config "$TESTDIR/../data/auspice_config1.json" \
  >   --output dataset.json &>/dev/null

  $ deep diff --ignore-order  "$TESTDIR/../data/dataset1.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']"
  {}

...same but with repeated --node-data options instead of a single multi-valued option

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --node-data "$TESTDIR/../data/div_node-data.json" \
  >   --node-data "$TESTDIR/../data/location_node-data.json" \
  >   --auspice-config "$TESTDIR/../data/auspice_config1.json" \
  >   --output dataset.json &>/dev/null

  $ deep diff --ignore-order  "$TESTDIR/../data/dataset1.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']"
  {}
