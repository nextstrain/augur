Setup

  $ source "$TESTDIR"/_setup.sh

Metadata with ambiguous days on all strains should error when grouping by week.

  $ cat >metadata.tsv <<~~
  > strain	date
  > SEQ1	2000-01-XX
  > SEQ2	2000-02-XX
  > SEQ3	2000-03-XX
  > SEQ4	2000-04-XX
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by week \
  >   --sequences-per-group 1 \
  >   --subsample-seed 0 \
  >   --output-metadata metadata-filtered.tsv \
  >   --output-log filtered_log.tsv
  4 strains were dropped during filtering
  	4 were dropped during grouping due to ambiguous day information
  	0 were dropped because of subsampling criteria
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  [2]
  $ cat filtered_log.tsv | grep "skip_group_by_with_ambiguous_day" | wc -l
  \s*4 (re)
  $ cat metadata-filtered.tsv
  strain	date

Metadata with ambiguous months on all strains should error when grouping by month.

  $ cat >metadata.tsv <<~~
  > strain	date
  > SEQ1	2000-XX-XX
  > SEQ2	2000-XX-XX
  > SEQ3	2000-XX-XX
  > SEQ4	2000-XX-XX
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by month \
  >   --sequences-per-group 1 \
  >   --subsample-seed 0 \
  >   --output-metadata metadata-filtered.tsv \
  >   --output-log filtered_log.tsv
  4 strains were dropped during filtering
  	4 were dropped during grouping due to ambiguous month information
  	0 were dropped because of subsampling criteria
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  [2]
  $ cat filtered_log.tsv | grep "skip_group_by_with_ambiguous_month" | wc -l
  \s*4 (re)
  $ cat metadata-filtered.tsv
  strain	date

Metadata with ambiguous years on all strains should error when grouping by year.

  $ cat >metadata.tsv <<~~
  > strain	date
  > SEQ1	XXXX-XX-XX
  > SEQ2	XXXX-XX-XX
  > SEQ3	XXXX-XX-XX
  > SEQ4	XXXX-XX-XX
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --group-by year \
  >   --sequences-per-group 1 \
  >   --subsample-seed 0 \
  >   --output-metadata metadata-filtered.tsv \
  >   --output-log filtered_log.tsv
  4 strains were dropped during filtering
  	4 were dropped during grouping due to ambiguous year information
  	0 were dropped because of subsampling criteria
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  [2]
  $ cat filtered_log.tsv | grep "skip_group_by_with_ambiguous_year" | wc -l
  \s*4 (re)
  $ cat metadata-filtered.tsv
  strain	date
