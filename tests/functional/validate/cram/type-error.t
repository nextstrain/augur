Setup

  $ source "$TESTDIR"/_setup.sh

Test wrong data type validation:

  $ cat >wrong-type.json <<~~
  > {
  >   "version": "v2",
  >   "tree": {
  >     "name": "test_tree",
  >     "node_attrs": "should_be_object_not_string"
  >   }
  > }
  > ~~

Type mismatch error:

  $ ${AUGUR} validate export-v2 wrong-type.json
  Validating schema of 'wrong-type.json'...
    .tree failed: {"name": "test_tree", "node_attrs": "should_be_o…} is invalid, see below.
      validation for arm 0: {"$ref": "#/$defs/tree"}
        .tree.node_attrs failed: Expected object but found str '"should_be_object_not_string"'
      validation for arm 1: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree failed: Expected array but found dict '{"name": "test_tree", "node_attrs": "should_be_o…}'
    top level failed: Missing required field 'meta'
  FATAL ERROR: Validation of 'wrong-type.json' failed.
  [2]
