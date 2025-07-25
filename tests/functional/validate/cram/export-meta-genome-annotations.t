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
    .meta.genome_annotations.invalid_gene "not_an_object" failed oneOf validation for [{"$ref": "#/$defs/startend"}, {"$ref": "#/$defs/segments"}]
      validation for arm 0: {"$ref": "#/$defs/startend"}
        .meta.genome_annotations.invalid_gene "not_an_object" failed type validation for "object"
      validation for arm 1: {"$ref": "#/$defs/segments"}
        .meta.genome_annotations.invalid_gene "not_an_object" failed type validation for "object"
    .meta.genome_annotations.invalid_gene "not_an_object" failed type validation for "object"
    .meta.genome_annotations {"invalid_gene": "not_an_object"} failed required validation for ["nuc"]
  FATAL ERROR: Validation of 'invalid-annotations.json' failed.
  [2]
