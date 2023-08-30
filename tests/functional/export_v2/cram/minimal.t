Setup

  $ source "$TESTDIR"/_setup.sh

Minimal export

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --node-data "$TESTDIR/../data/div_node-data.json" \
  >   --output minimal.json
  WARNING: You didn't provide information on who is maintaining this analysis.
  
  Validating produced JSON
  Validating schema of 'minimal.json'...
  Validating that the JSON is internally consistent...
  Validation of 'minimal.json' succeeded.
  

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py"  "$TESTDIR/../data/minimal.json" minimal.json \
  >   --exclude-paths "root['meta']['updated']"
  {}

Future test:
Run augur export _without_ any node-data JSONs when this can read divergence values from the newick file
and compare this to the tree using `div_node-data.json` - they should be identical.
See https://github.com/nextstrain/augur/issues/273 for more
