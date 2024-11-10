Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata file for testing.

  $ cat >metadata.tsv <<~~
  > strain	coverage	category
  > SEQ_1	0.94	A
  > SEQ_2	0.95	B
  > SEQ_3	0.96	C
  > SEQ_4		
  > ~~

Automatic inference works.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "coverage >= 0.95 & category == 'B'" \
  >  --output-strains filtered_strains.txt
  3 strains were dropped during filtering
  	3 were filtered out by the query: "coverage >= 0.95 & category == 'B'"
  1 strain passed all filters

Specifying coverage:float explicitly also works.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "coverage >= 0.95 & category == 'B'" \
  >  --query-columns coverage:float \
  >  --output-strains filtered_strains.txt
  3 strains were dropped during filtering
  	3 were filtered out by the query: "coverage >= 0.95 & category == 'B'"
  1 strain passed all filters

Specifying coverage:float category:str also works.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "coverage >= 0.95 & category == 'B'" \
  >  --query-columns coverage:float category:str \
  >  --output-strains filtered_strains.txt
  3 strains were dropped during filtering
  \t3 were filtered out by the query: "coverage >= 0.95 & category == 'B'" (esc)
  1 strain passed all filters

Specifying category:float does not work.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "coverage >= 0.95 & category == 'B'" \
  >  --query-columns category:float \
  >  --output-strains filtered_strains.txt
  ERROR: Failed to convert value in column 'category' to float. Unable to parse string "A" at position 0
  [2]

Specifying category:bool also does not work.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query "coverage >= 0.95 & category == 'B'" \
  >  --query-columns category:bool \
  >  --output-strains filtered_strains.txt
  ERROR: Failed to convert value in column 'category' to bool. Unable to convert 'A' to a boolean value.
  [2]
