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
    .tree failed: {"name": "root", "node_attrs": {"div": 0, "num_d…} is invalid, see below.
      .tree.children[…].node_attrs.invalid_nested_structure failed: Expected object but found bool 'true'
      .tree.children[…].node_attrs failed: {"invalid_nested_structure": true} is invalid, see below.
        .tree.children[…].node_attrs failed: Missing required field 'div'
        .tree.children[…].node_attrs failed: Missing required field 'num_date'
      .tree.children[…].node_attrs.num_date failed: Expected object but found float '2023.2'
      .tree.children[0].node_attrs.num_date failed: Expected object but found float '2023.1'
      .tree.node_attrs.num_date failed: Expected object but found float '2023.0'
      .tree failed: Expected array but found dict '{"name": "root", "node_attrs": {"div": 0, "num_d\…}'
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
    .meta.genome_annotations.invalid_gene failed: "not_an_object" is invalid, see below.
      .meta.genome_annotations.invalid_gene failed: Expected object but found str '"not_an_object"'
      .meta.genome_annotations.invalid_gene failed: Expected object but found str '"not_an_object"'
    .meta.genome_annotations.invalid_gene failed: Expected object but found str '"not_an_object"'
    .meta.genome_annotations failed: Missing required field 'nuc'
    .tree failed: {"name": "root", "node_attrs": {"div": 0, "num_d…} is invalid, see below.
      .tree.node_attrs.num_date failed: Expected object but found float '2023.0'
      .tree failed: Expected array but found dict '{"name": "root", "node_attrs": {"div": 0, "num_d\…}'
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
    .tree failed: {"name": "root", "node_attrs": {"div": 0, "num_d…} is invalid, see below.
      .tree.node_attrs.num_date failed: Expected object but found float '2023.0'
      .tree failed: Expected array but found dict '{"name": "root", "node_attrs": {"div": 0, "num_d\…}'
  FATAL ERROR: Validation of 'mixed-refs.json' failed.
  [2]
