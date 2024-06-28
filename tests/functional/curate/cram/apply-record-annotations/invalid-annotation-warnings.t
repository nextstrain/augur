Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Using invalid annotations results in warning messages.
Valid annotations are still applied.

  $ cat >annotations.tsv <<~~
  > record_1	field_1	annotation_1	extra_field
  > record_1	field_1	annotation_1
  > record_1	field_1
  > record_1
  > ~~

  $ echo '{"accession": "record_1", "field_1": "value_1"}' \
  >   |  ${AUGUR} curate apply-record-annotations \
  >       --annotations annotations.tsv
  WARNING: Could not decode annotation line record_1\tfield_1\tannotation_1\textra_field (esc)
  WARNING: Could not decode annotation line record_1\tfield_1 (esc)
  WARNING: Could not decode annotation line record_1
  {"accession": "record_1", "field_1": "annotation_1"}
