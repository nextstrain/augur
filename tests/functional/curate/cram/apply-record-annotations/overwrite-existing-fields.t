Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Test that annotations overwrite existing fields.

  $ cat >annotations.tsv <<~~
  > record_1	field_1	annotation_1
  > record_1	field_2	annotation_2
  > ~~

  $ echo '{"accession": "record_1", "field_1": "value_1", "field_2": "value_2"}' \
  >   |  ${AUGUR} curate apply-record-annotations \
  >       --annotations annotations.tsv
  {"accession": "record_1", "field_1": "annotation_1", "field_2": "annotation_2"}
