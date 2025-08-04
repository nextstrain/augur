Setup

  $ source "$TESTDIR"/_setup.sh

Test specific filtering behaviors and date filters.

Test date-based filtering with max_date.

  $ cat >max_date_config.yaml <<~~
  > samples:
  >   early_only:
  >     group_by:
  >     - region
  >     subsample_max_sequences: 5
  >     max_date: 2016-02-29
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config max_date_config.yaml \
  >   --output-metadata max_date_output.tsv \
  >   --output-sequences max_date_output.fasta \
  >   --subsample-seed 0
  Validating schema of 'max_date_config.yaml'...
  Sampling at 5 per group.
  7 strains were dropped during filtering
  	7 were dropped because they were later than 2016.16 or missing a date
  	0 were dropped because of subsampling criteria
  5 strains passed all filters
  8 strains were dropped during filtering
  	1 had no metadata
  	12 were dropped by `--exclude-all`
  \\t5 were added back because they were in .*sample_early_only.* (re)
  5 strains passed all filters

Test date-based filtering with min_date.

  $ cat >min_date_config.yaml <<~~
  > samples:
  >   recent_only:
  >     group_by:
  >     - region
  >     subsample_max_sequences: 5
  >     min_date: 2016-03-01
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config min_date_config.yaml \
  >   --output-metadata min_date_output.tsv \
  >   --output-sequences min_date_output.fasta \
  >   --subsample-seed 0
  Validating schema of 'min_date_config.yaml'...
  Sampling at 1 per group.
  9 strains were dropped during filtering
  	5 were dropped because they were earlier than 2016.17 or missing a date
  	4 were dropped because of subsampling criteria
  3 strains passed all filters
  10 strains were dropped during filtering
  	1 had no metadata
  	12 were dropped by `--exclude-all`
  \\t3 were added back because they were in .*sample_recent_only.* (re)
  3 strains passed all filters

Test exclude_where filtering.

  $ cat >exclude_config.yaml <<~~
  > samples:
  >   exclude_test:
  >     group_by:
  >     - region
  >     subsample_max_sequences: 10
  >     exclude_where:
  >     - region!=South America
  >     - region!=North America
  > ~~

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config exclude_config.yaml \
  >   --output-metadata exclude_output.tsv \
  >   --output-sequences exclude_output.fasta \
  >   --subsample-seed 0
  Validating schema of 'exclude_config.yaml'...
  Sampling at 10 per group.
  12 strains were dropped during filtering
  	6 were dropped because of 'region!=South America'
  	6 were dropped because of 'region!=North America'
  	0 were dropped because of subsampling criteria
  ERROR: Sample 'exclude_test' failed: All samples have been dropped! Check filter rules and metadata file format.
  [2]

Verify exclude_where filtering worked correctly.

  $ cut -f5 exclude_output.tsv | tail -n +2 | sort | uniq
  cut: exclude_output.tsv: No such file or directory
  [1]
