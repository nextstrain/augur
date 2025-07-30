Setup

  $ source "$TESTDIR"/_setup.sh

Test that lengthy paths are truncated but still readable.

  $ cat >nested-error.json <<~~
  > {
  >   "version": "v2",
  >   "tree": {
  >     "name": "root",
  >     "node_attrs": {},
  >     "children": [
  >       {
  >         "name": "child1",
  >         "node_attrs": {},
  >         "children": [
  >           {
  >             "name": "grandchild1",
  >             "node_attrs": {},
  >             "children": [
  >               {
  >                 "name": "great_grandchild",
  >                 "invalid_field": "this_should_not_be_here"
  >               }
  >             ]
  >           }
  >         ]
  >       }
  >     ]
  >   }
  > }
  > ~~

  $ ${AUGUR} validate export-v2 nested-error.json
  Validating schema of 'nested-error.json'...
    .tree failed: {"name": "root", "node_attrs": {}, "children": […} did not match one of the acceptable options below.
      Option 1: {"$ref": "#/$defs/tree"}
        .tree.children[…] failed: Unexpected property 'invalid_field'
        .tree.children[…] failed: Missing required property 'node_attrs'
        .tree.children[…].node_attrs failed: {} did not match any of the acceptable options below.
          Option 1: {"required": ["div"]}
            .tree.children[…].node_attrs failed: Missing required property 'div'
          Option 2: {"required": ["num_date"]}
            .tree.children[…].node_attrs failed: Missing required property 'num_date'
        .tree.children[0].node_attrs failed: {} did not match any of the acceptable options below.
          Option 1: {"required": ["div"]}
            .tree.children[0].node_attrs failed: Missing required property 'div'
          Option 2: {"required": ["num_date"]}
            .tree.children[0].node_attrs failed: Missing required property 'num_date'
        .tree.node_attrs failed: {} did not match any of the acceptable options below.
          Option 1: {"required": ["div"]}
            .tree.node_attrs failed: Missing required property 'div'
          Option 2: {"required": ["num_date"]}
            .tree.node_attrs failed: Missing required property 'num_date'
      Option 2: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree failed: Expected array but found dict '{"name": "root", "node_attrs": {}, "children": […}'
    top level failed: Missing required property 'meta'
  FATAL ERROR: Validation of 'nested-error.json' failed.
  [2]
