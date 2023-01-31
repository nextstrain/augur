Setup

  $ source "$TESTDIR"/_setup.sh

Testing combined metadata and FASTA inputs for the curate command.
Running the `passthru` subcommand since it does not do any data transformations.

Create FASTA file for testing.

  $ cat >sequences.fasta <<~~
  > >sequence_A
  > ATCG
  > >sequence_B
  > TCGA
  > >sequence_C
  > CGAT
  > ~~

Create metadata TSV file for testing.

  $ cat >metadata.tsv <<~~
  > strain	country	date
  > sequence_A	USA	2020-10-01
  > sequence_B	USA	2020-10-02
  > sequence_C	USA	2020-10-03
  > ~~

Test metadata input with extra FASTA input options without a FASTA file.
This is expected to fail with an error.

  $ ${AUGUR} curate passthru \
  > --metadata metadata.tsv \
  > --seq-id-column name \
  > --seq-field sequences
  ERROR: The --seq-id-column and --seq-field options should only be used when providing a FASTA file.
  [2]


Test metadata and FASTA inputs without required FASTA input options.
This is expected to fail with an error.

  $ ${AUGUR} curate passthru \
  > --metadata metadata.tsv \
  > --fasta sequences.fasta
  ERROR: The --seq-id-column and --seq-field options are required for a FASTA file input.
  [2]

Test metadata and FASTA inputs with required FASTA input options.

  $ ${AUGUR} curate passthru \
  > --metadata metadata.tsv \
  > --fasta sequences.fasta \
  > --seq-id-column strain \
  > --seq-field seq
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01", "seq": "ATCG"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02", "seq": "TCGA"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03", "seq": "CGAT"}

Create new metadata file with duplicate and extra metadata records.

  $ cp metadata.tsv metadata-with-duplicate-and-unmatched-records.tsv
  $ cat >>metadata-with-duplicate-and-unmatched-records.tsv <<~~
  > sequence_A	USA	2020-10-XX
  > extra_metadata_A	USA	2020-10-01
  > extra_metadata_B	USA	2020-10-02
  > ~~

Create new FASTA file with duplicate and extra sequence records.

  $ cp sequences.fasta sequences-with-duplicate-and-unmatched-records.fasta
  $ cat >>sequences-with-duplicate-and-unmatched-records.fasta <<~~
  > >sequence_A
  > NNNN
  > >extra_sequence_A
  > ATCG
  > >extra_sequence_B
  > TCGA
  > ~~

Test metadata and FASTA inputs with duplicate and extra records and default `ERROR_FIRST` reporting.
This is expected to fail with an error, so redirecting stdout since we don't care about the output.

  $ ${AUGUR} curate passthru \
  > --metadata metadata-with-duplicate-and-unmatched-records.tsv \
  > --fasta sequences-with-duplicate-and-unmatched-records.fasta \
  > --seq-id-column strain \
  > --seq-field seq 1> /dev/null
  ERROR: Encountered sequence record with duplicate id 'sequence_A'.
  [2]

Test metadata and FASTA inputs with duplicate and extra records with `ERROR_ALL` reporting.
This is expected to fail with an error, so redirecting stdout since we don't care about the output.

  $ ${AUGUR} curate passthru \
  > --metadata metadata-with-duplicate-and-unmatched-records.tsv \
  > --fasta sequences-with-duplicate-and-unmatched-records.fasta \
  > --seq-id-column strain \
  > --seq-field seq \
  > --unmatched-reporting error_all \
  > --duplicate-reporting error_all 1> /dev/null
  ERROR: Encountered the following error(s) when parsing metadata with sequences:
  The output may not match expectations because there were records with duplicate sequence ids.
  The following sequence ids were duplicated in .*metadata-with-duplicate-and-unmatched-records.* (re)
  'sequence_A'
  The following sequence ids were duplicated in .*sequences-with-duplicate-and-unmatched-records.* (re)
  'sequence_A'
  The output may be incomplete because there were unmatched records.
  The following metadata records did not have a matching sequence:
  'extra_metadata_A'
  'extra_metadata_B'
  The following sequence records did not have a matching metadata record:
  'extra_sequence_A'
  'extra_sequence_B'
  [2]

Test metadata and FASTA inputs with unmatched records, but ask to only warn on unmatched and duplicates.
This is expected run without error and only print a warning.
Notice the duplicate sequence "sequence_A" will always use the first sequence in the FASTA file because of pyfastx.

  $ ${AUGUR} curate passthru \
  > --metadata metadata-with-duplicate-and-unmatched-records.tsv \
  > --fasta sequences-with-duplicate-and-unmatched-records.fasta \
  > --seq-id-column strain \
  > --seq-field seq \
  > --unmatched-reporting warn \
  > --duplicate-reporting warn
  WARNING: Encountered sequence record with duplicate id 'sequence_A'.
  WARNING: Encountered metadata record with duplicate id 'sequence_A'.
  WARNING: Encountered metadata record 'extra_metadata_A' without a matching sequence.
  WARNING: Encountered metadata record 'extra_metadata_B' without a matching sequence.
  WARNING: The output may not match expectations because there were records with duplicate sequence ids.
  The following sequence ids were duplicated in .*metadata-with-duplicate-and-unmatched-records.* (re)
  'sequence_A'
  The following sequence ids were duplicated in .*sequences-with-duplicate-and-unmatched-records.* (re)
  'sequence_A'
  WARNING: The output may be incomplete because there were unmatched records.
  The following metadata records did not have a matching sequence:
  'extra_metadata_A'
  'extra_metadata_B'
  The following sequence records did not have a matching metadata record:
  'extra_sequence_A'
  'extra_sequence_B'
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01", "seq": "ATCG"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02", "seq": "TCGA"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03", "seq": "CGAT"}
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-XX", "seq": "ATCG"}

Test metadata and FASTA inputs with unmatched records in both, but ask to silent unmatched and duplicates.
Notice the duplicate sequence "sequence_A" will always use the first sequence in the FASTA file because of pyfastx.

  $ ${AUGUR} curate passthru \
  > --metadata metadata-with-duplicate-and-unmatched-records.tsv \
  > --fasta sequences-with-duplicate-and-unmatched-records.fasta \
  > --seq-id-column strain \
  > --seq-field seq \
  > --unmatched-reporting silent \
  > --duplicate-reporting silent
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-01", "seq": "ATCG"}
  {"strain": "sequence_B", "country": "USA", "date": "2020-10-02", "seq": "TCGA"}
  {"strain": "sequence_C", "country": "USA", "date": "2020-10-03", "seq": "CGAT"}
  {"strain": "sequence_A", "country": "USA", "date": "2020-10-XX", "seq": "ATCG"}

