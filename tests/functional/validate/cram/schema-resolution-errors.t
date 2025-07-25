Setup

  $ source "$TESTDIR"/_setup.sh

Test edge cases that might trigger schema resolution issues:

Test deeply nested structure:

  $ cat >deep-nested.json <<~~
  > {
  >   "version": "v2",
  >   "meta": {
  >     "updated": "2023-12-25",
  >     "panels": ["tree"],
  >     "genome_annotations": {
  >       "nuc": {
  >         "start": 1,
  >         "end": 1000,
  >         "strand": "+",
  >         "type": "source"
  >       }
  >     }
  >   },
  >   "tree": {
  >     "name": "root",
  >     "node_attrs": {"div": 0, "num_date": 2023.0},
  >     "children": [
  >       {
  >         "name": "child1",
  >         "node_attrs": {"div": 0.1, "num_date": 2023.1},
  >         "children": [
  >           {
  >             "name": "grandchild1",
  >             "node_attrs": {"div": 0.2, "num_date": 2023.2},
  >             "children": [
  >               {
  >                 "name": "great_grandchild1",
  >                 "node_attrs": {"invalid_nested_structure": true}
  >               }
  >             ]
  >           }
  >         ]
  >       }
  >     ]
  >   }
  > }
  > ~~

  $ ${AUGUR} validate export-v2 deep-nested.json
  Validating schema of 'deep-nested.json'...
    .tree {"name": "root", "node_attrs": {"div": 0, "num_d…} failed oneOf validation for [{"$ref": "#/$defs/tree"}, {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}]
      validation for arm 0: {"$ref": "#/$defs/tree"}
        .tree.children[…].node_attrs.invalid_nested_structure true failed type validation for "object"
        .tree.children[…].node_attrs {"invalid_nested_structure": true} failed anyOf validation for [{"required": ["div"]}, {"required": ["num_date"]}]
          validation for arm 0: {"required": ["div"]}
            .tree.children[…].node_attrs {"invalid_nested_structure": true} failed required validation for ["div"]
          validation for arm 1: {"required": ["num_date"]}
            .tree.children[…].node_attrs {"invalid_nested_structure": true} failed required validation for ["num_date"]
        .tree.children[…].node_attrs.num_date 2023.2 failed type validation for "object"
        .tree.children[0].node_attrs.num_date 2023.1 failed type validation for "object"
        .tree.node_attrs.num_date 2023.0 failed type validation for "object"
      validation for arm 1: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree {"name": "root", "node_attrs": {"div": 0, "num_d…} failed type validation for "array"
  FATAL ERROR: Validation of 'deep-nested.json' failed.
  [2]

Test with invalid genome annotations structure:

  $ cat >invalid-annotations.json <<~~
  > {
  >   "version": "v2",
  >   "meta": {
  >     "updated": "2023-12-25",
  >     "panels": ["tree"],
  >     "genome_annotations": {
  >       "invalid_gene": "not_an_object"
  >     }
  >   },
  >   "tree": {
  >     "name": "root",
  >     "node_attrs": {"div": 0, "num_date": 2023.0}
  >   }
  > }
  > ~~

  $ ${AUGUR} validate export-v2 invalid-annotations.json
  Validating schema of 'invalid-annotations.json'...
    .meta.genome_annotations.invalid_gene "not_an_object" failed oneOf validation for [{"$ref": "#/$defs/startend"}, {"$ref": "#/$defs/segments"}]
      validation for arm 0: {"$ref": "#/$defs/startend"}
        .meta.genome_annotations.invalid_gene "not_an_object" failed type validation for "object"
      validation for arm 1: {"$ref": "#/$defs/segments"}
        .meta.genome_annotations.invalid_gene "not_an_object" failed type validation for "object"
    .meta.genome_annotations.invalid_gene "not_an_object" failed type validation for "object"
    .meta.genome_annotations {"invalid_gene": "not_an_object"} failed required validation for ["nuc"]
    .tree {"name": "root", "node_attrs": {"div": 0, "num_d…} failed oneOf validation for [{"$ref": "#/$defs/tree"}, {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}]
      validation for arm 0: {"$ref": "#/$defs/tree"}
        .tree.node_attrs.num_date 2023.0 failed type validation for "object"
      validation for arm 1: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree {"name": "root", "node_attrs": {"div": 0, "num_d…} failed type validation for "array"
  FATAL ERROR: Validation of 'invalid-annotations.json' failed.
  [2]

Test with mixed valid and invalid references:

  $ cat >mixed-refs.json <<~~
  > {
  >   "version": "v2",
  >   "meta": {
  >     "updated": "2023-12-25",
  >     "panels": ["tree"],
  >     "maintainers": [
  >       {
  >         "name": "Valid Author",
  >         "url": "invalid-url-format"
  >       }
  >     ]
  >   },
  >   "tree": {
  >     "name": "root",
  >     "node_attrs": {"div": 0, "num_date": 2023.0}
  >   }
  > }
  > ~~

  $ ${AUGUR} validate export-v2 mixed-refs.json
  Validating schema of 'mixed-refs.json'...
    .tree {"name": "root", "node_attrs": {"div": 0, "num_d…} failed oneOf validation for [{"$ref": "#/$defs/tree"}, {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}]
      validation for arm 0: {"$ref": "#/$defs/tree"}
        .tree.node_attrs.num_date 2023.0 failed type validation for "object"
      validation for arm 1: {"type": "array", "minItems": 1, "items": {"$ref": "#/$defs/tree"}}
        .tree {"name": "root", "node_attrs": {"div": 0, "num_d…} failed type validation for "array"
  FATAL ERROR: Validation of 'mixed-refs.json' failed.
  [2]
