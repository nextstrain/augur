Integration tests for augur subsample.

  $ source "$TESTDIR"/_setup.sh

Create a simple subsampling config.

  $ cat >samples.yaml <<~~
  > samples:
  >   early:
  >     group_by:
  >     - region
  >     max_sequences: 1
  >     max_date: 2016-02-29
  >     exclude_where:
  >     - region!=South America
  >     subsample_seed: 0
  >   recent:
  >     group_by:
  >     - region
  >     max_sequences: 1
  >     min_date: 2016-03-01
  >     exclude_where:
  >     - region!=South America
  >     subsample_seed: 0
  > ~~

Run the subsample command.

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR/filter/data/metadata.tsv" \
  >   --sequences "$TESTDIR/filter/data/sequences.fasta" \
  >   --config samples.yaml \
  >   --output-metadata subsampled.tsv \
  >   --output-sequences subsampled.fasta
  Sampling at 1 per group.
  11 strains were dropped during filtering
  \t6 were dropped because of 'region!=South America' (esc)
  \t4 were dropped because they were later than 2016.16 or missing a date (esc)
  \t1 was dropped because of subsampling criteria (esc)
  1 strain passed all filters
  Sampling at 1 per group.
  11 strains were dropped during filtering
  \t6 were dropped because of 'region!=South America' (esc)
  \t3 were dropped because they were earlier than 2016.17 or missing a date (esc)
  \t2 were dropped because of subsampling criteria (esc)
  1 strain passed all filters
  11 strains were dropped during filtering
  \t1 had no metadata (esc)
  \t12 were dropped by `--exclude-all` (esc)
  \t1 was added back because it was in sample-early.txt (esc)
  \t1 was added back because it was in sample-recent.txt (esc)
  2 strains passed all filters

Check that two sequences remain and all come from South America.

  $ grep -c '^>' subsampled.fasta
  2
  $ tail -n +2 subsampled.tsv | wc -l | tr -d ' '
  2
  $ cut -f5 subsampled.tsv | tail -n +2 | sort | uniq
  South America
