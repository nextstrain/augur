Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Test that annotations for new fields are added to the record.

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
  {"accession": "record_1", "field_1": "value_1"}
  {"accession": "record_2", "field_1": "value_1", "new_field": "annotation_1"}


  $ cat records.ndjson \
  >   |  ${AUGUR} curate apply-record-annotations \
  >       --annotations annotations.tsv \
  >       --output-metadata metadata.tsv

  $ cat metadata.tsv
  accession\tfield_1 (esc)
  record_1\tvalue_1 (esc)
  record_2\tvalue_1 (esc)
