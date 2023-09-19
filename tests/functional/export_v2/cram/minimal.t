Setup

  $ source "$TESTDIR"/_setup.sh

Minimal export -- single input (tree) and single output (dataset JSON)

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --output minimal.json
  WARNING: You didn't provide information on who is maintaining this analysis.
  
  Validating produced JSON
  Validating schema of 'minimal.json'...
  Validating that the JSON is internally consistent...
  Validation of 'minimal.json' succeeded.
  

The above minimal.json takes divergence from the newick file. This converts newick divergences of (e.g.) '1' to `1.0`
because BioPython uses floats (which is perfectly reasonable). Ignore this type change in the JSON diff.
(Note that Auspice won't behave any differently)

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" --ignore-numeric-type-changes "$TESTDIR/../data/minimal.json" minimal.json \
  >   --exclude-paths "root['meta']['updated']"
  {}

Almost minimal export -- divergence is encoded via the node-data JSON typically produced by `augur refine`

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --node-data "$TESTDIR/../data/div_node-data.json" \
  >   --output almost-minimal.json
  WARNING: You didn't provide information on who is maintaining this analysis.
  
  Validating produced JSON
  Validating schema of 'almost-minimal.json'...
  Validating that the JSON is internally consistent...
  Validation of 'almost-minimal.json' succeeded.
  

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py"  "$TESTDIR/../data/minimal.json" almost-minimal.json \
  >   --exclude-paths "root['meta']['updated']"
  {}

Future test:
Run augur export _without_ any node-data JSONs when this can read divergence values from the newick file
and compare this to the tree using `div_node-data.json` - they should be identical.
See https://github.com/nextstrain/augur/issues/273 for more
