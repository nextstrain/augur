Setup

  $ source "$TESTDIR"/_setup.sh

Test with missing node_attrs and meta.

  $ cat >missing-required.json <<~~
  > {
  >   "version": "v2",
  >   "tree": {
  >     "name": "test_tree"
  >   }
  > }
  > ~~

  $ ${AUGUR} validate export-v2 missing-required.json
  Validating schema of 'missing-required.json'...
    .tree {"name": "test_tree"} did not match one of the acceptable options below.
      Option 1: {"$ref": "#/$defs/tree"}
        .tree {"name": "test_tree"} failed required validation for ["name", "node_attrs"]
      Option 2: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree {"name": "test_tree"} failed type validation for "array"
     {"version": "v2", "tree": {"name": "test_tree"}} failed required validation for ["version", "meta", "tree"]
  FATAL ERROR: Validation of 'missing-required.json' failed.
  [2]
