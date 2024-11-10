Setup

  $ source "$TESTDIR"/_setup.sh

Check output of min/max date filters.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --min-date 2015-01-01 \
  >  --max-date 2016-02-01 \
  >  --output-metadata filtered_metadata.tsv
  8 strains were dropped during filtering
  	1 was dropped because it was earlier than 2015.0 or missing a date
  	7 were dropped because they were later than 2016.09 or missing a date
  4 strains passed all filters
