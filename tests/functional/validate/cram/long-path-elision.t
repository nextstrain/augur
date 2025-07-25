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
    .tree {"name": "root", "node_attrs": {}, "children": […} failed oneOf validation for [{"$ref": "#/$defs/tree"}, {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}]
      validation for arm 0: {"$ref": "#/$defs/tree"}
        .tree.children[…] {"name": "great_grandchild", "invalid_field": "t…} failed additionalProperties validation for false
        .tree.children[…] {"name": "great_grandchild", "invalid_field": "t…} failed required validation for ["name", "node_attrs"]
        .tree.children[…].node_attrs {} failed anyOf validation for [{"required": ["div"]}, {"required": ["num_date"]}]
          validation for arm 0: {"required": ["div"]}
            .tree.children[…].node_attrs {} failed required validation for ["div"]
          validation for arm 1: {"required": ["num_date"]}
            .tree.children[…].node_attrs {} failed required validation for ["num_date"]
        .tree.children[0].node_attrs {} failed anyOf validation for [{"required": ["div"]}, {"required": ["num_date"]}]
          validation for arm 0: {"required": ["div"]}
            .tree.children[0].node_attrs {} failed required validation for ["div"]
          validation for arm 1: {"required": ["num_date"]}
            .tree.children[0].node_attrs {} failed required validation for ["num_date"]
        .tree.node_attrs {} failed anyOf validation for [{"required": ["div"]}, {"required": ["num_date"]}]
          validation for arm 0: {"required": ["div"]}
            .tree.node_attrs {} failed required validation for ["div"]
          validation for arm 1: {"required": ["num_date"]}
            .tree.node_attrs {} failed required validation for ["num_date"]
      validation for arm 1: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree {"name": "root", "node_attrs": {}, "children": […} failed type validation for "array"
     {"version": "v2", "tree": {"name": "root", "node…} failed required validation for ["version", "meta", "tree"]
  FATAL ERROR: Validation of 'nested-error.json' failed.
  [2]
