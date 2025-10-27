Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata TSV file with two date values that should be functionally equivalent.

  $ cat >metadata.tsv <<~~
  > strain	date
  > SEQ_1	
  > SEQ_2	XXXX-XX-XX
  > ~~

BUG: SEQ_2 passes for --min-date and --max-date.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --min-date 2025 \
  >  --output-strains filtered_strains.txt
  2 strains were dropped during filtering
  	2 were dropped because they were earlier than 2025.0 or missing a date
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  [2]

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --max-date 2025 \
  >  --output-strains filtered_strains.txt
  2 strains were dropped during filtering
  	2 were dropped because they were later than 2025.0 or missing a date
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  [2]
