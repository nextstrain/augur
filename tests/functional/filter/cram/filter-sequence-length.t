Setup

  $ source "$TESTDIR"/_setup.sh

Filter using --min-length and --max-length.

  $ ${AUGUR} filter \
  >  --sequence-index "$TESTDIR/../data/sequence_index.tsv" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --min-length 10500 \
  >  --max-length 10700 \
  >  --output-strains filtered_strains.txt
  7 strains were dropped during filtering
  	1 had no metadata
  	1 had no sequence data
  	2 were dropped because they were shorter than the minimum length of 10500bp when only counting standard nucleotide characters A, C, G, or T (case-insensitive)
  	3 were dropped because they were longer than the maximum length of 10700bp when only counting standard nucleotide characters A, C, G, or T (case-insensitive)
  6 strains passed all filters
