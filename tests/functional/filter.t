Integration tests for augur filter.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="../../bin/augur"

Filter with subsampling, requesting no more than 10 sequences.
With 10 groups to subsample from, this should produce one sequence per group.

  $ ${AUGUR} filter \
  >  --sequences filter/sequences.fasta \
  >  --sequence-index filter/sequences.fasta.idx \
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
  >  --sequence-index filter/sequences.fasta.idx \
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
  >  --sequence-index filter/sequences.fasta.idx \
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
  >  --sequence-index filter/sequences.fasta.idx \
  >  --metadata filter/metadata.tsv \
  >  --min-date 2012 \
  >  --group-by country year month \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 314159 \
  >  --output "$TMP/filtered.fasta" > /dev/null
  $ rm -f "$TMP/filtered.fasta"
