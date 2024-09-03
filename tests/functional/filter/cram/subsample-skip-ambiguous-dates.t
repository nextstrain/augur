Setup

  $ source "$TESTDIR"/_setup.sh

Try to subsample a maximum number of sequences by year and month, given metadata with ambiguous year and month values.
Strains with ambiguous years or months should be dropped and logged.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --group-by year month \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 0 \
  >  --output-strains filtered_strains.txt \
  >  --output-log filtered_log.tsv > /dev/null
  WARNING: Asked to provide at most 5 sequences, but there are 6 groups.
  Sampling probabilistically at 0.8333 sequences per group, meaning it is possible to have more than the requested maximum of 5 sequences after filtering.
  8 strains were dropped during filtering
  	1 was dropped during grouping due to ambiguous year information
  	1 was dropped during grouping due to ambiguous month information
  	6 were dropped because of subsampling criteria
  4 strains passed all filters
  $ grep "SG_018" filtered_log.tsv | cut -f 1-2
  SG_018\tskip_group_by_with_ambiguous_month (esc)
  $ grep "COL/FLR_00024/2015" filtered_log.tsv | cut -f 1-2
  COL/FLR_00024/2015\tskip_group_by_with_ambiguous_year (esc)

Group by 'year month week'. Using 'week' has some restrictions - 'year' should warn and 'month' should error.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --group-by year month week \
  >  --sequences-per-group 1 \
  >  --output-strains filtered_strains.txt > /dev/null
  WARNING: 'year' grouping will be ignored since 'week' includes ISO year.
  ERROR: 'month' and 'week' grouping cannot be used together.
  [2]

Group by 'week'. Check the number of strains that have been dropped due to ambiguous day.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --group-by week \
  >  --sequences-per-group 1 \
  >  --subsample-seed 0 \
  >  --output-strains filtered_strains.txt \
  >  --output-log filtered_log.tsv 2>/dev/null
  $ grep "skip_group_by_with_ambiguous_year" filtered_log.tsv | wc -l
  \s*1 (re)
  $ grep "skip_group_by_with_ambiguous_month" filtered_log.tsv | wc -l
  \s*1 (re)
  $ grep "skip_group_by_with_ambiguous_day" filtered_log.tsv | wc -l
  \s*3 (re)
