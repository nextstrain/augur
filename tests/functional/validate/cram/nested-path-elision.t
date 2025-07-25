Setup

  $ source "$TESTDIR"/_setup.sh

Test path elision in deeply nested JSON:

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

Path elision with children[…] format:

  $ ${AUGUR} validate export-v2 nested-error.json
  Validating schema of 'nested-error.json'...
    .tree {"name": "root", "node_attrs": {}, "children": […} is invalid, see below.
      validation for arm 0: {"$ref": "#/$defs/tree"}
        .tree.children[…] {"name": "great_grandchild", "invalid_field": "t…} failed: Additional properties are not allowed ('invalid_field' was unexpected)
        .tree.children[…] {"name": "great_grandchild", "invalid_field": "t…} failed: 'node_attrs' is a required property
        .tree.children[…].node_attrs {} failed: {} is not valid under any of the given schemas
          validation for arm 0: {"required": ["div"]}
            .tree.children[…].node_attrs {} failed: 'div' is a required property
          validation for arm 1: {"required": ["num_date"]}
            .tree.children[…].node_attrs {} failed: 'num_date' is a required property
        .tree.children[0].node_attrs {} failed: {} is not valid under any of the given schemas
          validation for arm 0: {"required": ["div"]}
            .tree.children[0].node_attrs {} failed: 'div' is a required property
          validation for arm 1: {"required": ["num_date"]}
            .tree.children[0].node_attrs {} failed: 'num_date' is a required property
        .tree.node_attrs {} failed: {} is not valid under any of the given schemas
          validation for arm 0: {"required": ["div"]}
            .tree.node_attrs {} failed: 'div' is a required property
          validation for arm 1: {"required": ["num_date"]}
            .tree.node_attrs {} failed: 'num_date' is a required property
      validation for arm 1: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree {"name": "root", "node_attrs": {}, "children": […} failed: {'name': 'root', 'node_attrs': {}, 'children': [{'name': 'child1', 'node_attrs': {}, 'children': [{'name': 'grandchild1', 'node_attrs': {}, 'children': [{'name': 'great_grandchild', 'invalid_field': 'this_should_not_be_here'}]}]}]} is not of type 'array'
     {"version": "v2", "tree": {"name": "root", "node…} failed: 'meta' is a required property
  FATAL ERROR: Validation of 'nested-error.json' failed.
  [2]
