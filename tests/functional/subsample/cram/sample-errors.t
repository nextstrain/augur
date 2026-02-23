Setup

  $ source "$TESTDIR"/_setup.sh

A failed sample should terminate the process early.

As long as we run with a single thread, the first sample (the bad one) will fail and
in our handling of that failure we can cancel subsequent samples which haven't yet
started

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
  >   --nthreads 1 \
  >   --seed 0
  Validating schema of 'samples.yaml'...
  [bad_sample] ERROR: You must specify a number of sequences per group or maximum sequences to subsample.
  ERROR: Sample 'bad_sample' failed, see error above.
  [2]
