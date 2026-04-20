Setup

  $ source "$TESTDIR"/_setup.sh

Test basic subsampling functionality with a simple config.
This recreates the `filter-sequence-length.t` test.

  $ cat >samples.yaml <<~~
  > samples:
  >   sample1:
  >     min_length: 329
  > ~~

  $ ${AUGUR} subsample \
  >   --seq-type aa \
  >   --metadata "$TESTDIR"/../../index/HA1.tsv \
  >   --sequences "$TESTDIR"/../../index/HA1.fasta \
  >   --config samples.yaml \
  >   --output-metadata subsampled.tsv \
  >   --output-sequences subsampled.fasta \
  >   --seed 0
  Validating schema of 'samples.yaml'...
  [sample1] running as the only filter call necessary
  [sample1] 3 strains were dropped during filtering
  [sample1] 	3 were dropped because they were shorter than the minimum length of 329 when only counting valid characters
  [sample1] 3 strains passed all filters

Check that three sequences remain in outputs.

  $ grep -c '^>' subsampled.fasta
  3

The same run, not specifying `--seq-type aa` will run as 'nuc', and thus
drop all sequences because their ATGC (i.e. nucleotide) lengths are small

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../index/HA1.tsv \
  >   --sequences "$TESTDIR"/../../index/HA1.fasta \
  >   --config samples.yaml \
  >   --output-metadata subsampled.tsv \
  >   --output-sequences subsampled.fasta \
  >   --seed 0
  Validating schema of 'samples.yaml'...
  [sample1] running as the only filter call necessary
  [sample1] 6 strains were dropped during filtering
  [sample1] 	6 were dropped because they were shorter than the minimum length of 329 when only counting valid characters
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  [2]


Proximal sampling is not yet supported

  $ cat >samples_prox.yaml <<~~
  > samples:
  >   sample1:
  >     min_length: 329
  >   prox:
  >     focal_sample: sample1
  > ~~

  $ ${AUGUR} subsample \
  >   --seq-type aa \
  >   --metadata "$TESTDIR"/../../index/HA1.tsv \
  >   --sequences "$TESTDIR"/../../index/HA1.fasta \
  >   --config samples_prox.yaml \
  >   --output-metadata subsampled.tsv \
  >   --output-sequences subsampled.fasta \
  >   --seed 0
  Validating schema of 'samples_prox.yaml'...
  ERROR: Proximal sampling for AA sequences is not yet supported.
  [2]
