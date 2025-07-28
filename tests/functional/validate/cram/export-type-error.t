Setup

  $ source "$TESTDIR"/_setup.sh

Test what happens when a string is given instead of an object.

  $ cat >wrong-type.json <<~~
  > {
  >   "version": "v2",
  >   "tree": {
  >     "name": "test_tree",
  >     "node_attrs": "should_be_object_not_string"
  >   }
  > }
  > ~~

  $ ${AUGUR} validate export-v2 wrong-type.json
  Validating schema of 'wrong-type.json'...
    .tree {"name": "test_tree", "node_attrs": "should_be_o…} did not match one of the acceptable options below.
      Option 1: {"$ref": "#/$defs/tree"}
        .tree.node_attrs "should_be_object_not_string" failed type validation for "object"
      Option 2: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree {"name": "test_tree", "node_attrs": "should_be_o…} failed type validation for "array"
     {"version": "v2", "tree": {"name": "test_tree", …} failed required validation for ["version", "meta", "tree"]
  FATAL ERROR: Validation of 'wrong-type.json' failed.
  [2]
