Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Filter with exclude query for two regions that comprise all but one strain.
This filter should leave a single record from Oceania.
Force include one South American record by country to get two total records.

Pandas engine
-------------

  $ ${AUGUR} filter \
  >  --metadata filter/data/metadata.tsv \
  >  --exclude-where "region=South America" "region=North America" "region=Southeast Asia" \
  >  --include-where "country=Ecuador" \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  $ wc -l "$TMP/filtered_strains.txt"
  \s*2 .* (re)
  $ rm -f "$TMP/filtered_strains.txt"

SQLite engine
-------------

  $ ${AUGUR} filter --engine sqlite \
  >  --metadata filter/data/metadata.tsv \
  >  --exclude-where "region=South America" "region=North America" "region=Southeast Asia" \
  >  --include-where "country=Ecuador" \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  $ wc -l "$TMP/filtered_strains.txt"
  \s*2 .* (re)
  $ rm -f "$TMP/filtered_strains.txt"
