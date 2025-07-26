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
    .tree failed: {"name": "test_tree"} did not match one of the acceptable options below.
      Option 1: {"$ref": "#/$defs/tree"}
        .tree failed: Missing required property 'node_attrs'
      Option 2: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree failed: Expected array but found dict '{"name": "test_tree"}'
    top level failed: Missing required property 'meta'
  FATAL ERROR: Validation of 'missing-required.json' failed.
  [2]
