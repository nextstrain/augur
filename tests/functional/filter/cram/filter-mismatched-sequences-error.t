Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Try to filter with sequences that don't match any of the metadata.
This should produce no results because the intersection of metadata and sequences is empty.

  $ echo -e ">something\nATCG" > "$TMP/dummy.fasta"
  $ ${AUGUR} filter \
  >  --sequences "$TMP/dummy.fasta" \
  >  --metadata filter/data/metadata.tsv \
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
  >  --metadata filter/data/metadata.tsv \
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
  >  --metadata filter/data/metadata.tsv \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  [1]
  $ wc -l "$TMP/filtered_strains.txt"
  \s*0 .* (re)
  $ rm -f "$TMP/filtered_strains.txt"
