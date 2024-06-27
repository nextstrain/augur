Setup

  $ source "$TESTDIR"/_setup.sh

Testing metadata outputs with internal quotes for the curate command.
Running the `passthru` subcommand since it does not do any data transformations.

Create NDJSON with internal quotes

  $ cat >records.ndjson <<~~
  > {"strain": "sequence_A", "submitting_lab": "SRC VB \"Vector\", Molecular Biology of Genomes"}
  > ~~

Test passthru with output to TSV.
This should not add any quotes around the field with internal quotes.

  $ cat records.ndjson \
  >   | ${AUGUR} curate passthru \
  >     --output-metadata output-metadata.tsv

  $ cat output-metadata.tsv
  strain\tsubmitting_lab (esc)
  sequence_A\tSRC VB "Vector", Molecular Biology of Genomes (esc)

Run the output TSV through augur curate passthru again.
The new output should still be identical to the first output.

  $ ${AUGUR} curate passthru \
  > --metadata output-metadata.tsv \
  > --output-metadata output-metadata-2.tsv

  $ diff -u output-metadata.tsv output-metadata-2.tsv
