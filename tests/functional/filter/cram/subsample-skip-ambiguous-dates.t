Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Try to subsample a maximum number of sequences by year and month, given metadata with ambiguous year and month values.
Strains with ambiguous years or months should be dropped and logged.

  $ ${AUGUR} filter \
  >  --metadata filter/data/metadata.tsv \
  >  --group-by year month \
  >  --subsample-max-sequences 5 \
  >  --output-strains "$TMP/filtered_strains.txt" \
  >  --output-log "$TMP/filtered_log.tsv" > /dev/null
  WARNING: Asked to provide at most 5 sequences, but there are 6 groups.
  $ grep "SG_018" "$TMP/filtered_log.tsv" | cut -f 1-2
  SG_018\tskip_group_by_with_ambiguous_month (esc)
  $ grep "COL/FLR_00024/2015" "$TMP/filtered_log.tsv" | cut -f 1-2
  COL/FLR_00024/2015\tskip_group_by_with_ambiguous_year (esc)
