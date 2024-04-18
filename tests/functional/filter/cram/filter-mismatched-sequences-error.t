Setup

  $ source "$TESTDIR"/_setup.sh

Try to filter with sequences that don't match any of the metadata.
This should produce no results because the intersection of metadata and sequences is empty.

  $ echo -e ">something\nATCG" > dummy.fasta
  $ ${AUGUR} filter \
  >  --sequences dummy.fasta \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --min-length 4 \
  >  --max-date 2020-01-30 \
  >  --output-strains filtered_strains.txt > /dev/null
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  [2]
  $ wc -l filtered_strains.txt
  wc: filtered_strains.txt: open: No such file or directory
  [1]

Repeat with sequence and strain outputs. We should get the same results.

  $ ${AUGUR} filter \
  >  --sequences dummy.fasta \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --max-date 2020-01-30 \
  >  --output-strains filtered_strains.txt \
  >  --output-sequences filtered.fasta > /dev/null
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  [2]
  $ wc -l filtered_strains.txt
  wc: filtered_strains.txt: open: No such file or directory
  [1]
  $ grep "^>" filtered.fasta | wc -l
  grep: filtered.fasta: No such file or directory
         0

Repeat without any sequence-based filters.
Since we expect metadata to be filtered by presence of strains in input sequences, this should produce no results because the intersection of metadata and sequences is empty.

  $ ${AUGUR} filter \
  >  --sequences dummy.fasta \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --output-strains filtered_strains.txt > /dev/null
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  [2]
  $ wc -l filtered_strains.txt
  wc: filtered_strains.txt: open: No such file or directory
  [1]
