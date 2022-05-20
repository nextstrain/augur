Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Filter with subsampling, requesting no more than 8 sequences.
With 8 groups to subsample from (after filtering), this should produce one sequence per group.

  $ ${AUGUR} filter \
  >  --sequences filter/data/sequences.fasta \
  >  --sequence-index filter/data/sequence_index.tsv \
  >  --metadata filter/data/metadata.tsv \
  >  --min-date 2012 \
  >  --group-by country year month \
  >  --subsample-max-sequences 8 \
  >  --subsample-seed 314159 \
  >  --no-probabilistic-sampling \
  >  --output "$TMP/filtered.fasta" > /dev/null
  $ grep ">" "$TMP/filtered.fasta" | wc -l
  \s*8 (re)
  $ rm -f "$TMP/filtered.fasta"
