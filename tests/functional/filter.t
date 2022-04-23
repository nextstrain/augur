Integration tests for augur filter.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="../../bin/augur"

Filter with exclude query for two regions that comprise all but one strain.
This filter should leave a single record from Oceania.
Force include one South American record by country to get two total records.

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --exclude-where "region=South America" "region=North America" "region=Southeast Asia" \
  >  --include-where "country=Ecuador" \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  $ wc -l "$TMP/filtered_strains.txt"
  \s*2 .* (re)
  $ rm -f "$TMP/filtered_strains.txt"

Filter with subsampling, requesting 1 sequence per group (for a group with 4 distinct values).

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --group-by region \
  >  --sequences-per-group 1 \
  >  --subsample-seed 314159 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  $ wc -l "$TMP/filtered_strains.txt"
  \s*4 .* (re)

By setting the subsample seed above, we should guarantee that we get the same "random" strains as another run with the same command.

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --group-by region \
  >  --sequences-per-group 1 \
  >  --subsample-seed 314159 \
  >  --output-strains "$TMP/filtered_strains_repeated.txt" > /dev/null

  $ diff -u <(sort "$TMP/filtered_strains.txt") <(sort "$TMP/filtered_strains_repeated.txt")
  $ rm -f "$TMP/filtered_strains.txt" "$TMP/filtered_strains_repeated.txt"

Filter with subsampling, requesting no more than 8 sequences.
With 8 groups to subsample from (after filtering), this should produce one sequence per group.

  $ ${AUGUR} filter \
  >  --sequences filter/sequences.fasta \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --min-date 2012 \
  >  --group-by country year month \
  >  --subsample-max-sequences 8 \
  >  --subsample-seed 314159 \
  >  --no-probabilistic-sampling \
  >  --output "$TMP/filtered.fasta" > /dev/null
  $ grep ">" "$TMP/filtered.fasta" | wc -l
  \s*8 (re)
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
  >  --output-strains "$TMP/filtered_strains_probabilistic.txt" > /dev/null
  WARNING: Asked to provide at most 5 sequences, but there are 8 groups.

Using the default probabilistic subsampling, should work the same as the previous case.

  $ ${AUGUR} filter \
  >  --sequences filter/sequences.fasta \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --min-date 2012 \
  >  --group-by country year month \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 314159 \
  >  --output-strains "$TMP/filtered_strains_default.txt" > /dev/null
  WARNING: Asked to provide at most 5 sequences, but there are 8 groups.

By setting the subsample seed above, we should get the same results for both runs.

  $ diff -u <(sort "$TMP/filtered_strains_probabilistic.txt") <(sort "$TMP/filtered_strains_default.txt")
  $ rm -f "$TMP/filtered_strains_probabilistic.txt" "$TMP/filtered_strains_default.txt"

Check output of probabilistic sampling.

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --group-by region year month \
  >  --subsample-max-sequences 3 \
  >  --probabilistic-sampling \
  >  --subsample-seed 314159 \
  >  --output-metadata "$TMP/filtered_metadata.tsv"
  WARNING: Asked to provide at most 3 sequences, but there are 8 groups.
  Sampling probabilistically at 0.3633 sequences per group, meaning it is possible to have more than the requested maximum of 3 sequences after filtering.
  10 strains were dropped during filtering
  \t1 were dropped during grouping due to ambiguous year information (esc)
  \t1 were dropped during grouping due to ambiguous month information (esc)
  \t8 of these were dropped because of subsampling criteria, using seed 314159 (esc)
  2 strains passed all filters

Ensure probabilistic sampling is not used when unnecessary.

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --group-by region year month \
  >  --subsample-max-sequences 10 \
  >  --probabilistic-sampling \
  >  --subsample-seed 314159 \
  >  --output-metadata "$TMP/filtered_metadata.tsv"
  Sampling at 10 per group.
  2 strains were dropped during filtering
  \t1 were dropped during grouping due to ambiguous year information (esc)
  \t1 were dropped during grouping due to ambiguous month information (esc)
  \t0 of these were dropped because of subsampling criteria, using seed 314159 (esc)
  10 strains passed all filters

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

Repeat this filter without a sequence index.
We should get the same outputs without building a sequence index on the fly, because the exclude-all flag tells us we only want to force-include strains and skip all other filters.

  $ ${AUGUR} filter \
  >  --sequences filter/sequences.fasta \
  >  --metadata filter/metadata.tsv \
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
  >  --sequences filter/sequences.fasta \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --exclude "$TMP/filtered_strains.brazil.txt" "$TMP/filtered_strains.colombia.txt" \
  >  --output "$TMP/filtered.fasta" > /dev/null
  $ grep "^>" "$TMP/filtered.fasta" | wc -l
  \s*7 (re)
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
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  [1]
  $ wc -l "$TMP/filtered_strains.txt"
  \s*0 .* (re)
  $ grep "^>" "$TMP/filtered.fasta" | wc -l
  \s*0 (re)
  $ rm -f "$TMP/filtered_strains.txt"
  $ rm -f "$TMP/filtered.fasta"

Repeat without any sequence-based filters.
Since we expect metadata to be filtered by presence of strains in input sequences, this should produce no results because the intersection of metadata and sequences is empty.

  $ ${AUGUR} filter \
  >  --sequences "$TMP/dummy.fasta" \
  >  --metadata filter/metadata.tsv \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  [1]
  $ wc -l "$TMP/filtered_strains.txt"
  \s*0 .* (re)
  $ rm -f "$TMP/filtered_strains.txt"

Filter TB strains from VCF and save as a list of filtered strains.

  $ ${AUGUR} filter \
  >  --sequences filter/tb.vcf.gz \
  >  --metadata filter/tb_metadata.tsv \
  >  --min-date 2012 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  $ wc -l "$TMP/filtered_strains.txt"
  \s*3 .* (re)
  $ rm -f "$TMP/filtered_strains.txt"

Confirm that filtering omits strains without metadata or sequences.
The input sequences are missing one strain that is in the metadata.
The metadata are missing one strain that has a sequence.
The list of strains to include has one strain with no metadata/sequence and one strain with information that would have been filtered by country.
The sequence index has one strain with invalid nucleotides.
The query initially filters 3 strains from Colombia, one of which is added back by the include.

  $ ${AUGUR} filter \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --query "country != 'Colombia'" \
  >  --non-nucleotide \
  >  --exclude-ambiguous-dates-by year \
  >  --include filter/include.txt \
  >  --output-strains "$TMP/filtered_strains.txt" \
  >  --output-log "$TMP/filtered_log.tsv"
  5 strains were dropped during filtering
  \t1 had no metadata (esc)
  \t1 had no sequence data (esc)
  \t3 of these were filtered out by the query: "country != 'Colombia'" (esc)
  \t1 of these were dropped because they had non-nucleotide characters (esc)
  \t1 strains were added back because they were in filter/include.txt (esc)
  8 strains passed all filters

  $ diff -u <(sort -k 1,1 filter/filtered_log.tsv) <(sort -k 1,1 "$TMP/filtered_log.tsv")
  $ rm -f "$TMP/filtered_strains.txt"

Subsample one strain per year with priorities.
There are two years (2015 and 2016) represented in the metadata.
The two highest priority strains are in these two years.

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --group-by year \
  >  --priority filter/priorities.tsv \
  >  --sequences-per-group 1 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null

  $ diff -u <(sort -k 2,2rn -k 1,1 filter/priorities.tsv | head -n 2 | cut -f 1) <(sort -k 1,1 "$TMP/filtered_strains.txt")
  $ rm -f "$TMP/filtered_strains.txt"

Try to subsample a maximum number of sequences by year and month, given metadata with ambiguous year and month values.
Strains with ambiguous years or months should be dropped and logged.

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --group-by year month \
  >  --subsample-max-sequences 5 \
  >  --output-strains "$TMP/filtered_strains.txt" \
  >  --output-log "$TMP/filtered_log.tsv" > /dev/null
  WARNING: Asked to provide at most 5 sequences, but there are 6 groups.
  $ grep "SG_018" "$TMP/filtered_log.tsv" | cut -f 1-2
  SG_018\tskip_group_by_with_ambiguous_month (esc)
  $ grep "COL/FLR_00024/2015" "$TMP/filtered_log.tsv" | cut -f 1-2
  COL/FLR_00024/2015\tskip_group_by_with_ambiguous_year (esc)

Try to group data without any grouping arguments.
This should fail with a helpful error message.

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --group-by year month \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  ERROR: You must specify a number of sequences per group or maximum sequences to subsample.
  [1]

Filter out a sequence with invalid nucleotides.

  $ ${AUGUR} filter \
  >  --sequence-index filter/sequence_index.tsv \
  >  --metadata filter/metadata.tsv \
  >  --non-nucleotide \
  >  --output-metadata "$TMP/filtered_metadata.tsv" > /dev/null
  $ wc -l "$TMP/filtered_metadata.tsv"
  \s*11 .* (re)

Try a comma-delimited file.

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.csv \
  >  --exclude-where "region=South America" "region=North America" "region=Southeast Asia" \
  >  --include-where "country=Ecuador" \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  $ wc -l "$TMP/filtered_strains.txt"
  \s*2 .* (re)
  $ rm -f "$TMP/filtered_strains.txt"

Error on duplicates in metadata within same chunk.

  $ cat >$TMP/metadata-duplicates.tsv <<~~
  > strain	date
  > a	2010-10-10
  > a	2010-10-10
  > b	2010-10-10
  > c	2010-10-10
  > d	2010-10-10
  > ~~
  $ ${AUGUR} filter \
  >   --metadata $TMP/metadata-duplicates.tsv \
  >   --group-by year \
  >   --sequences-per-group 2 \
  >   --subsample-seed 0 \
  >   --metadata-chunk-size 10 \
  >   --output-metadata $TMP/metadata-filtered.tsv > /dev/null
  ERROR: Duplicate found in .* (re)
  [2]
  $ cat $TMP/metadata-filtered.tsv
  cat: .*: No such file or directory (re)
  [1]

Error on duplicates in metadata in separate chunks.

  $ ${AUGUR} filter \
  >   --metadata $TMP/metadata-duplicates.tsv \
  >   --group-by year \
  >   --sequences-per-group 2 \
  >   --subsample-seed 0 \
  >   --metadata-chunk-size 1 \
  >   --output-metadata $TMP/metadata-filtered.tsv > /dev/null
  ERROR: Duplicate found in .* (re)
  [2]
  $ cat $TMP/metadata-filtered.tsv
  cat: .*: No such file or directory (re)
  [1]
