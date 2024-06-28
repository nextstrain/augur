Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Test annotations that use a custom id field.

  $ cat >annotations.tsv <<~~
  > record_1	field_1	annotation_1
  > ~~

  $ echo '{"record_id": "record_1", "field_1": "value_1"}' \
  >   |  ${AUGUR} curate apply-record-annotations \
  >       --id-field record_id \
  >       --annotations annotations.tsv
  {"record_id": "record_1", "field_1": "annotation_1"}
