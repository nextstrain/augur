Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata file

  $ cat >metadata.tsv <<~~
  > strain	region
  > SEQ_1	A
  > SEQ_2	B
  > SEQ_3	C
  > ~~

Scenario 1: Run command with one --exclude-where flag and multiple values

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >  --exclude-where "region=A" "region=B" \
  >  --output-strains filtered_strains.txt 2>/dev/null

Both exclusions are applied.

  $ cat filtered_strains.txt
  SEQ_3


Scenario 2: Run command with two --exclude-where flags

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >  --exclude-where "region=A" \
  >  --exclude-where "region=B" \
  >  --output-strains filtered_strains.txt 2>/dev/null

Both exclusions are applied.

  $ cat filtered_strains.txt
  SEQ_3
