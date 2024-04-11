Setup

  $ AUGUR="${AUGUR:-$TESTDIR/../../../bin/augur}"

Write data files.

  $ cat >metadata.tsv <<~~
  > strain	date	region
  > SEQ1	2021-01-01	A
  > SEQ2	2021-01-02	A
  > SEQ3	2021-01-01	B
  > SEQ4	2021-01-02	B
  > SEQ5	2021-02-02	B
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
  > size: 3
  > samples:
  >   focal:
  >     query: region=='A'
  >     weight: 1
  >   context:
  >     query: region=='B'
  >     weight: 2
  > ~~

Exclude:

  $ cat >exclude.txt <<~~
  > SEQ4
  > ~~

Apply subsampling.

FIXME: Use better regex to ignore temp file paths but still show other relevant info.

  $ ${AUGUR} subsample \
  >  --metadata metadata.tsv \
  >  --sequences sequences.fasta \
  >  --config config.yaml \
  >  --exclude exclude.txt \
  >  --output-metadata subsampled-metadata.tsv \
  >  --output-sequences subsampled-sequences.fasta \
  >  --random-seed 0 >/dev/null 2>&1

  $ cat subsampled-metadata.tsv
  strain	date	region
  SEQ1	2021-01-01	A
  SEQ3	2021-01-01	B
  SEQ5	2021-02-02	B
