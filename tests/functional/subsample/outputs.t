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
  >     filter: >-
  >       --query "region=='A'"
  >       --subsample-max-sequences 1
  >   context:
  >     filter: >-
  >       --query "region=='B'"
  >       --subsample-max-sequences 2
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
  >  --subsample-seed 0 \
  >  --tmpdir .
  
  Augur filter for intermediate sample 'focal' (no dependencies)
  RUNNING augur filter --query "region=='A'" --subsample-max-sequences 1 --metadata metadata.tsv --sequences sequences.fasta --output-strains ./focal.samples.txt --subsample-seed 0
  
  Augur filter for intermediate sample 'context' (no dependencies)
  RUNNING augur filter --query "region=='B'" --subsample-max-sequences 2 --metadata metadata.tsv --sequences sequences.fasta --output-strains ./context.samples.txt --subsample-seed 0
  
  Augur filter for intermediate sample 'output' depends on focal, context
  RUNNING augur filter --exclude-all --include ./focal.samples.txt ./context.samples.txt --metadata metadata.tsv --sequences sequences.fasta --output-metadata subsampled-metadata.tsv --output-sequences subsampled-sequences.fasta --subsample-seed 0

  $ cat subsampled-metadata.tsv
  strain	date	region
  SEQ1	2021-01-01	A
  SEQ3\t2021-01-01\tB (esc)
  SEQ4\t2021-01-02\tB (esc)

  $ cat subsampled-sequences.fasta
  >SEQ1
  AAA
  >SEQ3
  TTT
  >SEQ4
  GGG
