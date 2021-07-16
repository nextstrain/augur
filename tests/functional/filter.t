Integration tests for augur filter.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="../../bin/augur"

Filter with exclude query for two regions that comprise all but one strain.
This filter should leave a single record from Oceania.
Force include one South American record by country to get two total records.

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --exclude-where "region=South America" "region=North America" \
  >  --include-where "country=Ecuador" \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  $ wc -l "$TMP/filtered_strains.txt"
  \s*2 .* (re)
  $ rm -f "$TMP/filtered_strains.txt"

Filter with subsampling, requesting 1 sequence per group (for a group with 3 distinct values).

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --group-by region \
  >  --sequences-per-group 1 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  $ wc -l "$TMP/filtered_strains.txt"
  \s*3 .* (re)
  $ rm -f "$TMP/filtered_strains.txt"

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

Filter with subsampling where no more than 5 sequences are requested and no groups are specified.
This generates a dummy category and subsamples from there. With no-probabilistic-sampling we expect exactly 5 sequences.

  $ ${AUGUR} filter \
  >  --sequences filter/sequences.fasta \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --min-date 2012 \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 314159 \
  >  --no-probabilistic-sampling \
  >  --output "$TMP/filtered.fasta" > /dev/null
  $ grep ">" "$TMP/filtered.fasta" | wc -l
  \s*5 (re)
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
  ERROR: Asked to provide at most 5 sequences, but there are 8 groups.
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

Output should include the 8 sequences matching the filters and a header line.

  $ wc -l "$TMP/filtered_metadata.tsv"
  \s*9 .* (re)
  $ rm -f "$TMP/filtered_metadata.tsv"

Filter using only metadata and save results as a list of filtered strains.

  $ ${AUGUR} filter \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --min-date 2012 \
  >  --min-length 10500 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null

Output should include only the 8 sequences matching the filters (without a header line).

  $ wc -l "$TMP/filtered_strains.txt"
  \s*8 .* (re)
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
  \s*6 (re)
  $ rm -f "$TMP/filtered.fasta"

Try to filter with sequences that don't match any of the metadata.
This should produce no results because the intersection of metadata and sequences is empty.

  $ echo -e ">something\nATCG" > "$TMP/dummy.fasta"
  $ ${AUGUR} filter \
  >  --sequences "$TMP/dummy.fasta" \
  >  --metadata filter/metadata.tsv \
  >  --min-length 4 \
  >  --max-date 2020-01-30 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  [1]
  $ wc -l "$TMP/filtered_strains.txt"
  \s*0 .* (re)
  $ rm -f "$TMP/filtered_strains.txt"

Repeat with sequence and strain outputs. We should get the same results.

  $ ${AUGUR} filter \
  >  --sequences "$TMP/dummy.fasta" \
  >  --metadata filter/metadata.tsv \
  >  --max-date 2020-01-30 \
  >  --output-strains "$TMP/filtered_strains.txt" \
  >  --output-sequences "$TMP/filtered.fasta" > /dev/null
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  [1]
  $ wc -l "$TMP/filtered_strains.txt"
  \s*0 .* (re)
  $ grep "^>" "$TMP/filtered.fasta" | wc -l
  \s*0 (re)
  $ rm -f "$TMP/filtered_strains.txt"
  $ rm -f "$TMP/filtered.fasta"

Filter TB strains from VCF and save as a list of filtered strains.

  $ ${AUGUR} filter \
  >  --sequences filter/tb.vcf.gz \
  >  --metadata filter/tb_metadata.tsv \
  >  --min-date 2012 \
  >  --min-length 10500 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  $ wc -l "$TMP/filtered_strains.txt"
  \s*3 .* (re)
  $ rm -f "$TMP/filtered_strains.txt"

Confirm that filtering omits strains without metadata or sequences.
The input sequences are missing one strain that is in the metadata.
The metadata are missing one strain that has a sequence.
The list of strains to include has one strain with no metadata/sequence and one strain with information that would have been filtered by country.
The query initially filters 3 strains from Colombia, one of which is added back by the include.

  $ echo "NotReal" > "$TMP/include.txt"
  $ echo "COL/FLR_00008/2015" >> "$TMP/include.txt"
  $ ${AUGUR} filter \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --query "country != 'Colombia'" \
  >  --include "$TMP/include.txt" \
  >  --output-strains "$TMP/filtered_strains.txt"
  4 strains were dropped during filtering
  \t1 had no sequence data (esc)
  \t1 had no metadata (esc)
  \t3 of these were filtered out by the query: (esc)
  \t\t"country != 'Colombia'" (esc)
   (esc)
  \t1 strains were added back because they were requested by include files (esc)
  \t1 strains from include files were not added because they lacked sequence or metadata (esc)
  8 strains passed all filters

  $ rm -f "$TMP/filtered_strains.txt"
