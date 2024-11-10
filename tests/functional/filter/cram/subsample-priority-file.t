Setup

  $ source "$TESTDIR"/_setup.sh

Subsample one strain per year with priorities.
There are two years (2015 and 2016) represented in the metadata.
The two highest priority strains are in these two years.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --group-by year \
  >  --priority "$TESTDIR/../data/priorities.tsv" \
  >  --sequences-per-group 1 \
  >  --output-strains filtered_strains.txt 2>/dev/null

  $ diff -u <(sort -k 2,2rn -k 1,1 "$TESTDIR/../data/priorities.tsv" | head -n 2 | cut -f 1) <(sort -k 1,1 filtered_strains.txt)
