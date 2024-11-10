Setup

  $ source "$TESTDIR"/_setup.sh

Filter with exclude query for two regions that comprise all but one strain.
This filter should leave a single record from Oceania.
Force include one South American record by country to get two total records.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --exclude-where "region=South America" "region=North America" "region=Southeast Asia" \
  >  --include-where "country=Ecuador" \
  >  --output-strains filtered_strains.txt 2>/dev/null
  $ wc -l filtered_strains.txt
  \s*2 .* (re)
