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
    .tree {"name": "root", "node_attrs": {}, "children": […} did not match one of the acceptable options below.
      Option 1: {"$ref": "#/$defs/tree"}
        .tree.children[…] {"name": "great_grandchild", "invalid_field": "t…} failed: Additional properties are not allowed ('invalid_field' was unexpected)
        .tree.children[…] {"name": "great_grandchild", "invalid_field": "t…} failed: 'node_attrs' is a required property
        .tree.children[…].node_attrs {} did not match any of the acceptable options below.
          Option 1: {"required": ["div"]}
            .tree.children[…].node_attrs {} failed: 'div' is a required property
          Option 2: {"required": ["num_date"]}
            .tree.children[…].node_attrs {} failed: 'num_date' is a required property
        .tree.children[0].node_attrs {} did not match any of the acceptable options below.
          Option 1: {"required": ["div"]}
            .tree.children[0].node_attrs {} failed: 'div' is a required property
          Option 2: {"required": ["num_date"]}
            .tree.children[0].node_attrs {} failed: 'num_date' is a required property
        .tree.node_attrs {} did not match any of the acceptable options below.
          Option 1: {"required": ["div"]}
            .tree.node_attrs {} failed: 'div' is a required property
          Option 2: {"required": ["num_date"]}
            .tree.node_attrs {} failed: 'num_date' is a required property
      Option 2: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree {"name": "root", "node_attrs": {}, "children": […} failed: {'name': 'root', 'node_attrs': {}, 'children': [{'name': 'child1', 'node_attrs': {}, 'children': [{'name': 'grandchild1', 'node_attrs': {}, 'children': [{'name': 'great_grandchild', 'invalid_field': 'this_should_not_be_here'}]}]}]} is not of type 'array'
     {"version": "v2", "tree": {"name": "root", "node…} failed: 'meta' is a required property
  FATAL ERROR: Validation of 'nested-error.json' failed.
  [2]
