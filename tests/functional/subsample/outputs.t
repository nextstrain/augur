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

Apply subsampling.

  $ ${AUGUR} subsample \
  >  --metadata metadata.tsv \
  >  --sequences sequences.fasta \
  >  --config config.yaml \
  >  --output-metadata subsampled-metadata.tsv \
  >  --output-sequences subsampled-sequences.fasta \
  >  --random-seed 0
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  4 strains were dropped during filtering
  	3 were filtered out by the query: "region=='A'"
  	1 was dropped because of subsampling criteria
  1 strain passed all filters
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  3 strains were dropped during filtering
  	2 were filtered out by the query: "region=='B'"
  	1 was dropped because of subsampling criteria
  2 strains passed all filters
  2 strains were dropped during filtering
  	5 were dropped by `--exclude-all`
  .*1 was added back because it was in .*focal.samples.* (re)
  .*2 were added back because they were in .*context.samples.* (re)
  3 strains passed all filters
  Sampling for 'focal' (no dependencies)
  	query: region=='A'
  	max_sequences: 1
  
  augur filter .* (re)
  
  Sampling for 'context' (no dependencies)
  	query: region=='B'
  	max_sequences: 2
  
  augur filter .* (re)
  
  Sampling for 'output' (depends on focal, context)
  	exclude_all: True
  .*include:.* (re)
  
  augur filter .* (re)
  

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
