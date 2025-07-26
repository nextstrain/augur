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
    .tree failed: {"name": "root", "node_attrs": {}, "children": […} is invalid, see below.
      validation for arm 0: {"$ref": "#/$defs/tree"}
        .tree.children[…] failed: Unexpected property 'invalid_field'
        .tree.children[…] failed: Missing required field 'node_attrs'
        .tree.children[…].node_attrs failed: {} is invalid, see below.
          validation for arm 0: {"required": ["div"]}
            .tree.children[…].node_attrs failed: Missing required field 'div'
          validation for arm 1: {"required": ["num_date"]}
            .tree.children[…].node_attrs failed: Missing required field 'num_date'
        .tree.children[0].node_attrs failed: {} is invalid, see below.
          validation for arm 0: {"required": ["div"]}
            .tree.children[0].node_attrs failed: Missing required field 'div'
          validation for arm 1: {"required": ["num_date"]}
            .tree.children[0].node_attrs failed: Missing required field 'num_date'
        .tree.node_attrs failed: {} is invalid, see below.
          validation for arm 0: {"required": ["div"]}
            .tree.node_attrs failed: Missing required field 'div'
          validation for arm 1: {"required": ["num_date"]}
            .tree.node_attrs failed: Missing required field 'num_date'
      validation for arm 1: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree failed: Expected array but found dict '{"name": "root", "node_attrs": {}, "children": […}'
    top level failed: Missing required field 'meta'
  FATAL ERROR: Validation of 'nested-error.json' failed.
  [2]
