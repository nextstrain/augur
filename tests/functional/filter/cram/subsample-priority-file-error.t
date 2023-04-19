Setup

  $ source "$TESTDIR"/_setup.sh

Try running with a priority file that does not exist.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --priority does-not-exist.tsv \
  >  --subsample-max-sequences 5 \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: missing or malformed priority scores file does-not-exist.tsv
  [2]

Create a priority file that is not tab-delimited.

  $ cat >priorities.csv <<~~
  > SEQ_1,5
  > SEQ_2,6
  > SEQ_3,8
  > SEQ_4,-100
  > ~~

Try running with the above file.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --priority priorities.csv \
  >  --subsample-max-sequences 5 \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: missing or malformed priority scores file priorities.csv
  [2]

Create a priority file that does not have integers.

  $ cat >priorities.csv <<~~
  > SEQ_1	5
  > SEQ_2	6
  > SEQ_3	8
  > SEQ_4	a
  > ~~

Try running with the above file.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --priority priorities.csv \
  >  --subsample-max-sequences 5 \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: missing or malformed priority scores file priorities.csv
  [2]
