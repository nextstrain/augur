Setup

  $ source "$TESTDIR"/_setup.sh

  $ cat >metadata.tsv <<~~
  > strain	date
  > SEQ_1	
  > SEQ_2	2020
  > SEQ_3	2020-02-25
  > ~~

If we group by year month and some records don't have that information in
their date fields, we should skip those records from the group output and
track which records were skipped for which reasons.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --group-by year \
  >  --sequences-per-group 1 \
  >  --subsample-seed 0 \
  >  --output-log log.txt \
  >  --output-strains filtered_strains.txt > /dev/null
  2 strains were dropped during filtering
  	1 was dropped during grouping due to ambiguous year information
  	1 was dropped because of subsampling criteria
  1 strain passed all filters
  $ cat log.txt
  strain\tfilter\tkwargs (esc)
  SEQ_1\tskip_group_by_with_ambiguous_year\t"[[""date_column"", ""date""]]" (esc)
  SEQ_3\tsubsampling\t (esc)
  $ sort filtered_strains.txt
  SEQ_2

Similarly, if we group by month, we should skip records that don't have
month information in their date fields.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --group-by month \
  >  --sequences-per-group 1 \
  >  --subsample-seed 0 \
  >  --output-log log.txt \
  >  --output-strains filtered_strains.txt > /dev/null
  2 strains were dropped during filtering
  	1 was dropped during grouping due to ambiguous year information
  	1 was dropped during grouping due to ambiguous month information
  	0 were dropped because of subsampling criteria
  1 strain passed all filters
  $ cat log.txt
  strain\tfilter\tkwargs (esc)
  SEQ_1\tskip_group_by_with_ambiguous_year\t"[[""date_column"", ""date""]]" (esc)
  SEQ_2\tskip_group_by_with_ambiguous_month\t"[[""date_column"", ""date""]]" (esc)
  $ sort filtered_strains.txt
  SEQ_3
