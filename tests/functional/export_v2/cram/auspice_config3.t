Setup

  $ source "$TESTDIR"/_setup.sh

# auspice_config3.json is the same as auspice_config2.json but with an extra key which the schema does not allow.
# Running without --skip-validation should result in an error
# Message printed: "Validation of "$TESTDIR/../data/auspice_config3.json" failed."

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --node-data "$TESTDIR/../data/div_node-data.json" "$TESTDIR/../data/location_node-data.json" \
  >   --auspice-config "$TESTDIR/../data/auspice_config3.json" \
  >   --output dataset.json &>/dev/null
  [2]

# Skipping validation gives us the same results as `auspice_config2.json`

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --node-data "$TESTDIR/../data/div_node-data.json" "$TESTDIR/../data/location_node-data.json" \
  >   --auspice-config "$TESTDIR/../data/auspice_config3.json" \
  >   --output dataset.json \
  >   --skip-validation &>/dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py"  "$TESTDIR/../data/dataset2.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']"
  {}
