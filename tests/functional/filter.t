Integration tests for augur filter.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="../../bin/augur"

Filter with subsampling, requesting no more than 10 sequences.
With 10 groups to subsample from, this should produce one sequence per group.

  $ ${AUGUR} filter \
  >  --sequences filter/sequences.fasta \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --min-date 2012 \
  >  --group-by country year month \
  >  --subsample-max-sequences 10 \
  >  --subsample-seed 314159 \
  >  --no-probabilistic-sampling \
  >  --output "$TMP/filtered.fasta" > /dev/null
  $ grep ">" "$TMP/filtered.fasta" | wc -l
  \s*10 (re)
  $ rm -f "$TMP/filtered.fasta"

Try to filter with subsampling when there are more available groups than requested sequences.
This should fail, as probabilistic sampling is explicitly disabled.

  $ ${AUGUR} filter \
  >  --sequences filter/sequences.fasta \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --min-date 2012 \
  >  --group-by country year month \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 314159 \
  >  --no-probabilistic-sampling \
  >  --output "$TMP/filtered.fasta"
  ERROR: Asked to provide at most 5 sequences, but there are 10 groups.
  [1]
  $ rm -f "$TMP/filtered.fasta"

Explicitly use probabilistic subsampling to handle the case when there are more available groups than requested sequences.

  $ ${AUGUR} filter \
  >  --sequences filter/sequences.fasta \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --min-date 2012 \
  >  --group-by country year month \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 314159 \
  >  --probabilistic-sampling \
  >  --output "$TMP/filtered.fasta" > /dev/null
  $ rm -f "$TMP/filtered.fasta"

Using the default probabilistic subsampling, should work the same as the previous case.

  $ ${AUGUR} filter \
  >  --sequences filter/sequences.fasta \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --min-date 2012 \
  >  --group-by country year month \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 314159 \
  >  --output "$TMP/filtered.fasta" > /dev/null
  $ rm -f "$TMP/filtered.fasta"

Filter using only metadata without sequence input or output and save results as filtered metadata.

  $ ${AUGUR} filter \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --min-date 2012 \
  >  --min-length 10500 \
  >  --output-metadata "$TMP/filtered_metadata.tsv" > /dev/null

Output should include the 9 sequences matching the filters and a header line.

  $ wc -l "$TMP/filtered_metadata.tsv"
  \s*10 .* (re)
  $ rm -f "$TMP/filtered_metadata.tsv"

Filter using only metadata and save results as a list of filtered strains.

  $ ${AUGUR} filter \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --min-date 2012 \
  >  --min-length 10500 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null

Output should include only the 9 sequences matching the filters (without a header line).

  $ wc -l "$TMP/filtered_strains.txt"
  \s*9 .* (re)
  $ rm -f "$TMP/filtered_strains.txt"

Filter using only metadata without a sequence index.
This should work because the requested filters don't rely on sequence information.

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --min-date 2012 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  $ rm -f "$TMP/filtered_strains.txt"

Try to filter using only metadata without a sequence index.
This should fail because the requested filters rely on sequence information.

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --min-length 10000 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  ERROR: You need to provide a sequence index or sequences to filter on sequence-specific information.
  [1]

Try to filter with sequence outputs and no sequence inputs.
This should fail.

  $ ${AUGUR} filter \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --min-length 10000 \
  >  --output "$TMP/filtered.fasta" > /dev/null
  ERROR: You need to provide sequences to output sequences.
  [1]

Try to filter without any outputs.

  $ ${AUGUR} filter \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --min-length 10000 > /dev/null
  ERROR: You need to select at least one output.
  [1]

Filter into two separate sets and then select sequences from the union of those sets.
First, select strains from Brazil (there should be 1).

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --query "country == 'Brazil'" \
  >  --output-strains "$TMP/filtered_strains.brazil.txt" > /dev/null
  $ wc -l "$TMP/filtered_strains.brazil.txt"
  \s*1 .* (re)

Then, select strains from Colombia (there should be 3).

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --query "country == 'Colombia'" \
  >  --output-strains "$TMP/filtered_strains.colombia.txt" > /dev/null
  $ wc -l "$TMP/filtered_strains.colombia.txt"
  \s*3 .* (re)

Finally, exclude all sequences except those from the two sets of strains (there should be 4).

  $ ${AUGUR} filter \
  >  --sequences filter/sequences.fasta \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --exclude-all \
  >  --include "$TMP/filtered_strains.brazil.txt" "$TMP/filtered_strains.colombia.txt" \
  >  --output "$TMP/filtered.fasta" > /dev/null
  $ grep "^>" "$TMP/filtered.fasta" | wc -l
  \s*4 (re)
  $ rm -f "$TMP/filtered.fasta"

Alternately, exclude only the sequences from Brazil and Colombia (12 - 4 strains).

  $ ${AUGUR} filter \
  >  --sequences filter/sequences.fasta \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --exclude "$TMP/filtered_strains.brazil.txt" "$TMP/filtered_strains.colombia.txt" \
  >  --output "$TMP/filtered.fasta" > /dev/null
  $ grep "^>" "$TMP/filtered.fasta" | wc -l
  \s*8 (re)
  $ rm -f "$TMP/filtered.fasta"
