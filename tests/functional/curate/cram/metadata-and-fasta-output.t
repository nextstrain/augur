Setup

  $ source "$TESTDIR"/_setup.sh

Testing combined metadata and FASTA output for the curate command.
Running the `passthru` subcommand since it does not do any data transformations.

Create NDJSON file for testing.

  $ cat >records.ndjson <<~~
  > {"strain": "sequence_A", "country": "USA", "date": "2020-10-01", "sequence": "AAAA"}
  > {"strain": "sequence_T", "country": "USA", "date": "2020-10-02", "sequence": "TTTT"}
  > {"strain": "sequence_C", "country": "USA", "date": "2020-10-03", "sequence": "CCCC"}
  > ~~

Test metadata output with extra FASTA output options.
This is expected to fail immediately with an error.
  $ cat records.ndjson \
  >   | ${AUGUR} curate passthru \
  >     --output-metadata metadata.tsv \
  >     --output-id-field strain \
  >     --output-seq-field sequence
  ERROR: The --output-id-field and --output-seq-field options should only be used when requesting a FASTA output.
  [2]

Test metadata and FASTA outputs without requried FASTA output options.
This is expected to fail immediately with an error.
  $ cat records.ndjson \
  >   | ${AUGUR} curate passthru \
  >     --output-metadata metadata.tsv \
  >     --output-fasta sequences.fasta
  ERROR: The --output-id-field and --output-seq-field options are required for a FASTA output.
  [2]

Test metadata and FASTA outputs

  $ cat records.ndjson \
  >   | ${AUGUR} curate passthru \
  >     --output-metadata metadata.tsv \
  >     --output-fasta sequences.fasta \
  >     --output-id-field strain \
  >     --output-seq-field sequence
  $ cat metadata.tsv
  strain\tcountry\tdate (esc)
  sequence_A\tUSA\t2020-10-01 (esc)
  sequence_T\tUSA\t2020-10-02 (esc)
  sequence_C\tUSA\t2020-10-03 (esc)
  $ cat sequences.fasta
  >sequence_A (esc)
  AAAA (esc)
  >sequence_T (esc)
  TTTT (esc)
  >sequence_C (esc)
  CCCC

Test FASTA output without metadata output.

  $ cat records.ndjson \
  >   | ${AUGUR} curate passthru \
  >     --output-fasta sequences.fasta \
  >     --output-id-field strain \
  >     --output-seq-field sequence
  $ cat sequences.fasta
  >sequence_A (esc)
  AAAA (esc)
  >sequence_T (esc)
  TTTT (esc)
  >sequence_C (esc)
  CCCC

Test FASTA output with bad output id field.
This is expected to fail with an error.

  $ cat records.ndjson \
  >   | ${AUGUR} curate passthru \
  >     --output-fasta sequences.fasta \
  >     --output-id-field bogus_id \
  >     --output-seq-field sequence
  ERROR: Provided sequence identifier field 'bogus_id' does not exist.
  [2]

Test FASTA output with bad output sequence field.
This is expected to fail with an error.

  $ cat records.ndjson \
  >   | ${AUGUR} curate passthru \
  >     --output-fasta sequences.fasta \
  >     --output-id-field strain \
  >     --output-seq-field bogus_sequence
  ERROR: Provided sequence field 'bogus_sequence' does not exist.
  [2]
