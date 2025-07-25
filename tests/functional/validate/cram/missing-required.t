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
    .tree {"name": "test_tree"} failed oneOf validation for [{"$ref": "#/$defs/tree"}, {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}]
      validation for arm 0: {"$ref": "#/$defs/tree"}
        .tree {"name": "test_tree"} failed required validation for ["name", "node_attrs"]
      validation for arm 1: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree {"name": "test_tree"} failed type validation for "array"
     {"version": "v2", "tree": {"name": "test_tree"}} failed required validation for ["version", "meta", "tree"]
  FATAL ERROR: Validation of 'missing-required.json' failed.
  [2]
