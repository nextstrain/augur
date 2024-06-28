Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Using a custom id field that does not exist in the record results in an error.

  $ cat >annotations.tsv <<~~
  > record_1	field_1	annotation_1
  > ~~

  $ echo '{"record_id": "record_1", "field_1": "value_1"}' \
  >   |  ${AUGUR} curate apply-record-annotations \
  >       --id-field bad_id \
  >       --annotations annotations.tsv
  ERROR: ID field 'bad_id' does not exist in record
  [2]
