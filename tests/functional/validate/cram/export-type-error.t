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
    .tree failed: {"name": "test_tree", "node_attrs": "should_be_o…} did not match one of the acceptable options below.
      Option 1: {"$ref": "#/$defs/tree"}
        .tree.node_attrs failed: Expected object but found str '"should_be_object_not_string"'
      Option 2: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree failed: Expected array but found dict '{"name": "test_tree", "node_attrs": "should_be_o…}'
    top level failed: Missing required property 'meta'
  FATAL ERROR: Validation of 'wrong-type.json' failed.
  [2]
