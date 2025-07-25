Setup

  $ source "$TESTDIR"/_setup.sh

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
  >     "node_attrs": {"div": 0, "num_date": {"value": 2023.0} }
  >   }
  > }
  > ~~

  $ ${AUGUR} validate export-v2 invalid-annotations.json
  Validating schema of 'invalid-annotations.json'...
    .meta.genome_annotations.invalid_gene "not_an_object" did not match one of the acceptable options below.
      Option 1: {"$ref": "#/$defs/startend"}
        .meta.genome_annotations.invalid_gene "not_an_object" failed: 'not_an_object' is not of type 'object'
      Option 2: {"$ref": "#/$defs/segments"}
        .meta.genome_annotations.invalid_gene "not_an_object" failed: 'not_an_object' is not of type 'object'
    .meta.genome_annotations.invalid_gene "not_an_object" failed: 'not_an_object' is not of type 'object'
    .meta.genome_annotations {"invalid_gene": "not_an_object"} failed: 'nuc' is a required property
  FATAL ERROR: Validation of 'invalid-annotations.json' failed.
  [2]
