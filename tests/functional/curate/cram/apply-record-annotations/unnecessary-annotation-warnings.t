Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Using unnecessary annotations results in warning messages.
Valid annotations are still applied.

  $ cat >annotations.tsv <<~~
  > record_1	field_1	annotation_1
  > record_1	field_2	annotation_2
  > record_1	field_3	annotation_3
  > ~~

  $ echo '{"accession": "record_1", "field_1": "annotation_1", "field_2": "annotation_2", "field_3": "value_3"}' \
  >   |  ${AUGUR} curate apply-record-annotations \
  >       --annotations annotations.tsv
  WARNING: Unnecessary annotation for 'record_1' field 'field_1': value was already 'annotation_1'
  WARNING: Unnecessary annotation for 'record_1' field 'field_2': value was already 'annotation_2'
  {"accession": "record_1", "field_1": "annotation_1", "field_2": "annotation_2", "field_3": "annotation_3"}
