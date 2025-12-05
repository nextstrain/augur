Setup

  $ source "$TESTDIR"/_setup.sh

A failed sample should terminate the process early.

Note that good_sample_1 starts due to a race condition where the sample starts
before the previous sample's failure is detected.

  $ cat >samples.yaml <<~~
  > samples:
  >   bad_sample:
  >     group_by:
  >     - region
  >     max_sequences: 0
  >   good_sample_1:
  >     max_sequences: 1
  >   good_sample_2:
  >     max_sequences: 1
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config samples.yaml \
  >   --output-metadata subsampled.tsv \
  >   --output-sequences subsampled.fasta \
  >   --seed 0
  Validating schema of 'samples.yaml'...
  [bad_sample] ERROR: You must specify a number of sequences per group or maximum sequences to subsample.
  [good_sample_1] 11 strains were dropped during filtering
  [good_sample_1] 	11 were dropped because of subsampling criteria
  [good_sample_1] 1 strain passed all filters
  ERROR: Sample 'bad_sample' failed, see error above.
  [2]
