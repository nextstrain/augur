Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Test that annotations for new fields that are not in the record results in warning message.

  $ cat >annotations.tsv <<~~
  > record_2	new_field	annotation_1
  > ~~

  $ cat >records.ndjson <<~~
  > {"accession": "record_1", "field_1": "value_1"}
  > {"accession": "record_2", "field_1": "value_1"}
  > ~~

  $ cat records.ndjson \
  >   |  ${AUGUR} curate apply-record-annotations \
  >       --annotations annotations.tsv
  WARNING: Skipping annotation for field 'new_field' that does not exist in record
  {"accession": "record_1", "field_1": "value_1"}
  {"accession": "record_2", "field_1": "value_1"}
