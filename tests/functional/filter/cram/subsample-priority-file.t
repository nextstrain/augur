Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Subsample one strain per year with priorities.
There are two years (2015 and 2016) represented in the metadata.
The two highest priority strains are in these two years.

  $ ${AUGUR} filter \
  >  --metadata filter/data/metadata.tsv \
  >  --group-by year \
  >  --priority filter/data/priorities.tsv \
  >  --sequences-per-group 1 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null

  $ diff -u <(sort -k 2,2rn -k 1,1 filter/data/priorities.tsv | head -n 2 | cut -f 1) <(sort -k 1,1 "$TMP/filtered_strains.txt")
  $ rm -f "$TMP/filtered_strains.txt"
