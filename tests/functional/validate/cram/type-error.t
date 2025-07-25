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
    .tree {"name": "test_tree", "node_attrs": "should_be_o…} is invalid, see below.
      validation for arm 0: {"$ref": "#/$defs/tree"}
        .tree.node_attrs "should_be_object_not_string" failed: 'should_be_object_not_string' is not of type 'object'
      validation for arm 1: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree {"name": "test_tree", "node_attrs": "should_be_o…} failed: {'name': 'test_tree', 'node_attrs': 'should_be_object_not_string'} is not of type 'array'
     {"version": "v2", "tree": {"name": "test_tree", …} failed: 'meta' is a required property
  FATAL ERROR: Validation of 'wrong-type.json' failed.
  [2]
