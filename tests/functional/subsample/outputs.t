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
  >  --subsample-seed 0
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  4 strains were dropped during filtering
  \t3 were filtered out by the query: "region=='A'" (esc)
  \t1 was dropped because of subsampling criteria (esc)
  1 strain passed all filters
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  3 strains were dropped during filtering
  \t2 were filtered out by the query: "region=='B'" (esc)
  \t1 was dropped because of subsampling criteria (esc)
  2 strains passed all filters
  2 strains were dropped during filtering
  \t5 were dropped by `--exclude-all` (esc)
  \t1 was added back because it was in /var/folders/rk/s2p6nv1d13lcf0dynjkvxyl80000gp/T/cramtests-gum5cd6e/tmp/tmp6seok2rr/focal.samples.txt (esc)
  \t2 were added back because they were in /var/folders/rk/s2p6nv1d13lcf0dynjkvxyl80000gp/T/cramtests-gum5cd6e/tmp/tmp6seok2rr/context.samples.txt (esc)
  3 strains passed all filters
  
  Augur filter for intermediate sample 'focal' (no dependencies)
  RUNNING augur filter .* (re)
  
  Augur filter for intermediate sample 'context' (no dependencies)
  RUNNING augur filter .* (re)
  
  Augur filter for intermediate sample 'output' depends on focal, context
  RUNNING augur filter .* (re)

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
