Setup

  $ source "$TESTDIR"/_setup.sh

Run with a config file.

  $ cat > config.yaml <<~~
  > tree: "$TESTDIR/../data/tree.nwk"
  > node_data:
  >   - "$TESTDIR/../data/div_node-data.json"
  >   - "$TESTDIR/../data/location_node-data.json"
  > auspice_config: "$TESTDIR/../data/auspice_config1.json"
  > output: "dataset.json"
  > ~~

  $ ${AUGUR} export v2 --config config.yaml &>/dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset1.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']"
  {}
