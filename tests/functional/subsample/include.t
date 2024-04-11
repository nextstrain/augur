Setup

  $ AUGUR="${AUGUR:-$TESTDIR/../../../bin/augur}"

Write data files.

  $ cat >metadata.tsv <<~~
  > strain	date	region
  > SEQ1	2021-01-01	A
  > SEQ2	2021-01-02	A
  > SEQ3	2021-01-01	B
  > SEQ4	2021-01-02	B
  > SEQ5	2021-02-02	C
  > ~~

  $ cat >sequences.fasta <<~~
  > >SEQ1
  > AAA
  > >SEQ2
  > CCC
  > >SEQ3
  > TTT
  > >SEQ4
  > GGG
  > >SEQ5
  > AAC
  > ~~

Subsampling configuration:

  $ cat >config.yaml <<~~
  > samples:
  >   focal:
  >     query: region=='A'
  >     max_sequences: 1
  >   context:
  >     query: region=='B'
  >     max_sequences: 2
  > ~~

  $ cat >root.txt <<~~
  > SEQ5
  > ~~

Apply subsampling.

  $ ${AUGUR} subsample \
  >  --metadata metadata.tsv \
  >  --sequences sequences.fasta \
  >  --config config.yaml \
  >  --include root.txt \
  >  --output-metadata subsampled-metadata.tsv \
  >  --output-sequences subsampled-sequences.fasta \
  >  --random-seed 0 >/dev/null 2>&1

  $ cat subsampled-metadata.tsv
  strain	date	region
  SEQ1	2021-01-01	A
  SEQ3	2021-01-01	B
  SEQ4	2021-01-02	B
  SEQ5	2021-02-02	C
