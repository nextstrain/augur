Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Check output of min/max date filters.

Pandas engine
-------------

  $ ${AUGUR} filter \
  >  --metadata filter/data/metadata.tsv \
  >  --min-date 2015-01-01 \
  >  --max-date 2016-02-01 \
  >  --output-metadata "$TMP/filtered_metadata.tsv"
  8 strains were dropped during filtering
  \t1 of these were dropped because they were earlier than 2015.0 or missing a date (esc)
  \t7 of these were dropped because they were later than 2016.09 or missing a date (esc)
  4 strains passed all filters

SQLite engine
-------------

  $ ${AUGUR} filter --engine sqlite \
  >  --metadata filter/data/metadata.tsv \
  >  --min-date 2015-01-01 \
  >  --max-date 2016-02-01 \
  >  --output-metadata "$TMP/filtered_metadata.tsv"
  8 strains were dropped during filtering
  \t1 of these were dropped because they were earlier than 2015.0 or missing a date (esc)
  \t7 of these were dropped because they were later than 2016.09 or missing a date (esc)
  4 strains passed all filters
