Setup

  $ source "$TESTDIR"/_setup.sh

Test that node data attributes of different types are handled appropriately.

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --node-data "$TESTDIR/../data/node-data-types.json" \
  >   --output dataset.json &> /dev/null

Verify the expected node attrs from the node-data JSON are exported as filters
and colorings. Not diffing the dataset.json because order of node_attrs is not
guaranteed without an auspice_config.json

TODO: Switch to use jq once it's available in a well-defined test environment.
<https://github.com/nextstrain/augur/issues/1557>

  $ python3 -c 'import json, sys; print(sorted(json.load(sys.stdin)["meta"]["filters"]))' < dataset.json
  ['test_bool', 'test_float', 'test_int', 'test_str']

  $ python3 -c 'import json, sys; print(sorted(c["key"] for c in json.load(sys.stdin)["meta"]["colorings"]))' < dataset.json
  ['test_bool', 'test_float', 'test_int', 'test_str']
