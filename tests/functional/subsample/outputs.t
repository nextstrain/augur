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
  > samples:
  >   focal:
  >     filter:
  >       query: region=='A'
  >       subsample_max_sequences: 1
  >   context:
  >     filter:
  >       query: region=='B'
  >       subsample_max_sequences: 2
  > output:
  >   - focal
  >   - context
  > ~~

Apply subsampling.

  $ ${AUGUR} subsample \
  >  --config config.yaml \
  >  --metadata metadata.tsv \
  >  --sequences sequences.fasta \
  >  --output-metadata subsampled-metadata.tsv \
  >  --output-sequences subsampled-sequences.fasta \
  >  --subsample-seed 0
  4 strains were dropped during filtering
  	3 were filtered out by the query: "region=='A'"
  	1 was dropped because of subsampling criteria
  1 strain passed all filters
  3 strains were dropped during filtering
  	2 were filtered out by the query: "region=='B'"
  	1 was dropped because of subsampling criteria
  2 strains passed all filters
  2 strains were dropped during filtering
  	5 were dropped by `--exclude-all`
  .* (re)
  .* (re)
  3 strains passed all filters
  RUNNING augur filter with name 'focal' (no dependencies)
  	metadata: metadata.tsv
  .* (re)
  	query: region=='A'
  	subsample_max_sequences: 1
  
  RUNNING augur filter with name 'context' (no dependencies)
  	metadata: metadata.tsv
  .* (re)
  	query: region=='B'
  	subsample_max_sequences: 2
  
  RUNNING augur filter with name 'output' depends on focal, context
  	metadata: metadata.tsv
  	sequences: sequences.fasta
  	output_metadata: subsampled-metadata.tsv
  	output_sequences: subsampled-sequences.fasta
  	exclude_all: True
  .* (re)
  

  $ cat subsampled-metadata.tsv
  strain	date	region
  SEQ1	2021-01-01	A
  SEQ3	2021-01-01	B
  SEQ4	2021-01-02	B

  $ cat subsampled-sequences.fasta
  >SEQ1
  AAA
  >SEQ3
  TTT
  >SEQ4
  GGG
