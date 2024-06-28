Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../../bin/augur}"

Test that comments in the annotations file are ignored.

  $ cat >annotations.tsv <<~~
  > # This is a comment.
  > record_1	field_1	annotation_1 # This is also a comment.
  > ~~

  $ echo '{"accession": "record_1", "field_1": "value_1"}' \
  >   |  ${AUGUR} curate apply-record-annotations \
  >       --annotations annotations.tsv
  {"accession": "record_1", "field_1": "annotation_1"}
