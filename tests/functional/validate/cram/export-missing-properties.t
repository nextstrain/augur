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
        .tree {"name": "test_tree"} failed: 'node_attrs' is a required property
      Option 2: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree {"name": "test_tree"} failed: {'name': 'test_tree'} is not of type 'array'
     {"version": "v2", "tree": {"name": "test_tree"}} failed: 'meta' is a required property
  FATAL ERROR: Validation of 'missing-required.json' failed.
  [2]
