Setup

  $ source "$TESTDIR"/_setup.sh

Test duplicate node names in tree:

  $ cat >duplicate-nodes.json <<~~
  > {
  >   "version": "v2",
  >   "meta": {
  >     "updated": "2023-12-25",
  >     "panels": ["tree"]
  >   },
  >   "tree": {
  >     "name": "root",
  >     "node_attrs": {"div": 0, "num_date": 2023.0},
  >     "children": [
  >       {
  >         "name": "duplicate",
  >         "node_attrs": {"div": 0.1, "num_date": 2023.1}
  >       },
  >       {
  >         "name": "duplicate",
  >         "node_attrs": {"div": 0.2, "num_date": 2023.2}
  >       }
  >     ]
  >   }
  > }
  > ~~

  $ ${AUGUR} validate export-v2 duplicate-nodes.json
  Validating schema of 'duplicate-nodes.json'...
    .tree failed: {"name": "root", "node_attrs": {"div": 0, "num_d…} is invalid, see below.
      validation for arm 0: {"$ref": "#/$defs/tree"}
        .tree.children[0].node_attrs.num_date failed: Expected object but found float '2023.1'
        .tree.children[1].node_attrs.num_date failed: Expected object but found float '2023.2'
        .tree.node_attrs.num_date failed: Expected object but found float '2023.0'
      validation for arm 1: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree failed: Expected array but found dict '{"name": "root", "node_attrs": {"div": 0, "num_d…}'
  FATAL ERROR: Validation of 'duplicate-nodes.json' failed.
  [2]

Test export with no colorings:

  $ cat >no-colorings.json <<~~
  > {
  >   "version": "v2",
  >   "meta": {
  >     "updated": "2023-12-25",
  >     "panels": ["tree"]
  >   },
  >   "tree": {
  >     "name": "root",
  >     "node_attrs": {"div": 0, "num_date": 2023.0}
  >   }
  > }
  > ~~

  $ ${AUGUR} validate export-v2 no-colorings.json
  Validating schema of 'no-colorings.json'...
    .tree failed: {"name": "root", "node_attrs": {"div": 0, "num_d…} is invalid, see below.
      validation for arm 0: {"$ref": "#/$defs/tree"}
        .tree.node_attrs.num_date failed: Expected object but found float '2023.0'
      validation for arm 1: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree failed: Expected array but found dict '{"name": "root", "node_attrs": {"div": 0, "num_d…}'
  FATAL ERROR: Validation of 'no-colorings.json' failed.
  [2]

Test with missing colorings reference in meta:

  $ cat >missing-coloring-ref.json <<~~
  > {
  >   "version": "v2",
  >   "meta": {
  >     "updated": "2023-12-25",
  >     "panels": ["tree"],
  >     "colorings": [
  >       {
  >         "key": "region",
  >         "title": "Region"
  >       }
  >     ]
  >   },
  >   "tree": {
  >     "name": "root",
  >     "node_attrs": {"div": 0, "num_date": 2023.0},
  >     "children": [
  >       {
  >         "name": "child1",
  >         "node_attrs": {"div": 0.1, "num_date": 2023.1}
  >       }
  >     ]
  >   }
  > }
  > ~~

  $ ${AUGUR} validate export-v2 missing-coloring-ref.json
  Validating schema of 'missing-coloring-ref.json'...
    .meta.colorings[0] failed: Missing required field 'type'
    .tree failed: {"name": "root", "node_attrs": {"div": 0, "num_d…} is invalid, see below.
      validation for arm 0: {"$ref": "#/$defs/tree"}
        .tree.children[0].node_attrs.num_date failed: Expected object but found float '2023.1'
        .tree.node_attrs.num_date failed: Expected object but found float '2023.0'
      validation for arm 1: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree failed: Expected array but found dict '{"name": "root", "node_attrs": {"div": 0, "num_d…}'
  FATAL ERROR: Validation of 'missing-coloring-ref.json' failed.
  [2]
