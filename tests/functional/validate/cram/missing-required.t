Setup

  $ source "$TESTDIR"/_setup.sh

Testing that JSON schema validation errors produce human-readable messages
instead of technical schema validation output for issue #1426.

Create a minimal invalid JSON file missing required fields to test error message formatting.

  $ cat >missing-required.json <<~~
  > {
  >   "version": "v2",
  >   "tree": {
  >     "name": "test_tree"
  >   }
  > }
  > ~~

Missing required fields error:

  $ ${AUGUR} validate export-v2 missing-required.json
  Validating schema of 'missing-required.json'...
    .tree failed: {"name": "test_tree"} is invalid, see below.
      validation for arm 0: {"$ref": "#/$defs/tree"}
        .tree failed: Missing required field 'node_attrs'
      validation for arm 1: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree failed: Expected array but found dict '{"name": "test_tree"}'
    top level failed: Missing required field 'meta'
  FATAL ERROR: Validation of 'missing-required.json' failed.
  [2]
