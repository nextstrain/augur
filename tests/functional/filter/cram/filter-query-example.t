Setup

  $ source "$TESTDIR"/_setup.sh

Filter into two separate sets and then select sequences from the union of those sets.
First, select strains from Brazil (there should be 1).

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --query "country == 'Brazil'" \
  >  --output-strains filtered_strains.brazil.txt 2>/dev/null
  $ wc -l filtered_strains.brazil.txt
  \s*1 .* (re)

Then, select strains from Colombia (there should be 3).

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --query "country == 'Colombia'" \
  >  --output-strains filtered_strains.colombia.txt 2>/dev/null
  $ wc -l filtered_strains.colombia.txt
  \s*3 .* (re)

Finally, exclude all sequences except those from the two sets of strains (there should be 4).

  $ ${AUGUR} filter \
  >  --sequences "$TESTDIR/../data/sequences.fasta" \
  >  --sequence-index "$TESTDIR/../data/sequence_index.tsv" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --exclude-all \
  >  --include filtered_strains.brazil.txt filtered_strains.colombia.txt \
  >  --output filtered.fasta 2>/dev/null
  $ grep "^>" filtered.fasta | wc -l
  \s*4 (re)

Repeat this filter without a sequence index.
We should get the same outputs without building a sequence index on the fly, because the exclude-all flag tells us we only want to force-include strains and skip all other filters.

  $ ${AUGUR} filter \
  >  --sequences "$TESTDIR/../data/sequences.fasta" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --exclude-all \
  >  --include filtered_strains.brazil.txt filtered_strains.colombia.txt \
  >  --output filtered.fasta \
  >  --output-metadata filtered.tsv 2>/dev/null
  $ grep "^>" filtered.fasta | wc -l
  \s*4 (re)

Metadata should have the same number of records as the sequences plus a header.

  $ wc -l filtered.tsv
  \s*5 .* (re)

Alternately, exclude the sequences from Brazil and Colombia (N=4) and records without sequences (N=1) or metadata (N=1).

  $ ${AUGUR} filter \
  >  --sequences "$TESTDIR/../data/sequences.fasta" \
  >  --sequence-index "$TESTDIR/../data/sequence_index.tsv" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --exclude filtered_strains.brazil.txt filtered_strains.colombia.txt \
  >  --output filtered.fasta 2>/dev/null
  $ grep "^>" filtered.fasta | wc -l
  \s*7 (re)
