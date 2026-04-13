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

N (unknown) characters don't contribute to the length
Strain SG_018 has 10640 A/T/G/C + 19 Ns so check it's excluded from a 10,650 min length

  $ ${AUGUR} filter \
  >  --sequence-index "$TESTDIR/../data/sequence_index.tsv" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --min-length 10650 \
  >  --output-strains filtered_strains_2.txt \
  >  --output-log filtered_strains_2_log.txt &>/dev/null

  $ grep 'SG_018' filtered_strains_2_log.txt 
  SG_018\tfilter_by_min_length\t"[[""min_length"", 10650]]" (esc)


Other IUPAC characters don't contribute to the length either
Strain DOM/2016/BB_0059 has 9,408 A/T/G/C + 621 Ns + 6 others
Note: Ns are excluded form the length (see previous test) 

  $ ${AUGUR} filter \
  >  --sequence-index "$TESTDIR/../data/sequence_index.tsv" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --min-length 9410 \
  >  --output-strains filtered_strains_2.txt \
  >  --output-log filtered_strains_2_log.txt &>/dev/null

  $ grep 'DOM/2016/BB_0059' filtered_strains_2_log.txt 
  DOM/2016/BB_0059\tfilter_by_min_length\t"[[""min_length"", 9410]]" (esc)

