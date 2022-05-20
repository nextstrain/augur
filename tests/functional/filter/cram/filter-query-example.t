Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Filter into two separate sets and then select sequences from the union of those sets.
First, select strains from Brazil (there should be 1).

  $ ${AUGUR} filter \
  >  --metadata filter/data/metadata.tsv \
  >  --query "country == 'Brazil'" \
  >  --output-strains "$TMP/filtered_strains.brazil.txt" > /dev/null
  $ wc -l "$TMP/filtered_strains.brazil.txt"
  \s*1 .* (re)

Then, select strains from Colombia (there should be 3).

  $ ${AUGUR} filter \
  >  --metadata filter/data/metadata.tsv \
  >  --query "country == 'Colombia'" \
  >  --output-strains "$TMP/filtered_strains.colombia.txt" > /dev/null
  $ wc -l "$TMP/filtered_strains.colombia.txt"
  \s*3 .* (re)

Finally, exclude all sequences except those from the two sets of strains (there should be 4).

  $ ${AUGUR} filter \
  >  --sequences filter/data/sequences.fasta \
  >  --sequence-index filter/data/sequence_index.tsv \
  >  --metadata filter/data/metadata.tsv \
  >  --exclude-all \
  >  --include "$TMP/filtered_strains.brazil.txt" "$TMP/filtered_strains.colombia.txt" \
  >  --output "$TMP/filtered.fasta" > /dev/null
  $ grep "^>" "$TMP/filtered.fasta" | wc -l
  \s*4 (re)
  $ rm -f "$TMP/filtered.fasta"

Repeat this filter without a sequence index.
We should get the same outputs without building a sequence index on the fly, because the exclude-all flag tells us we only want to force-include strains and skip all other filters.

  $ ${AUGUR} filter \
  >  --sequences filter/data/sequences.fasta \
  >  --metadata filter/data/metadata.tsv \
  >  --exclude-all \
  >  --include "$TMP/filtered_strains.brazil.txt" "$TMP/filtered_strains.colombia.txt" \
  >  --output "$TMP/filtered.fasta" \
  >  --output-metadata "$TMP/filtered.tsv" > /dev/null
  $ grep "^>" "$TMP/filtered.fasta" | wc -l
  \s*4 (re)
  $ rm -f "$TMP/filtered.fasta"

Metadata should have the same number of records as the sequences plus a header.

  $ wc -l "$TMP/filtered.tsv"
  \s*5 .* (re)
  $ rm -f "$TMP/filtered.tsv"

Alternately, exclude the sequences from Brazil and Colombia (N=4) and records without sequences (N=1) or metadata (N=1).

  $ ${AUGUR} filter \
  >  --sequences filter/data/sequences.fasta \
  >  --sequence-index filter/data/sequence_index.tsv \
  >  --metadata filter/data/metadata.tsv \
  >  --exclude "$TMP/filtered_strains.brazil.txt" "$TMP/filtered_strains.colombia.txt" \
  >  --output "$TMP/filtered.fasta" > /dev/null
  $ grep "^>" "$TMP/filtered.fasta" | wc -l
  \s*7 (re)
  $ rm -f "$TMP/filtered.fasta"
