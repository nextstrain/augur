Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Explicitly use probabilistic subsampling to handle the case when there are more available groups than requested sequences.

  $ ${AUGUR} filter \
  >  --sequences filter/data/sequences.fasta \
  >  --sequence-index filter/data/sequence_index.tsv \
  >  --metadata filter/data/metadata.tsv \
  >  --min-date 2012 \
  >  --group-by country year month \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 314159 \
  >  --probabilistic-sampling \
  >  --output-strains "$TMP/filtered_strains_probabilistic.txt" > /dev/null
  WARNING: Asked to provide at most 5 sequences, but there are 8 groups.

Using the default probabilistic subsampling, should work the same as the previous case.

  $ ${AUGUR} filter \
  >  --sequences filter/data/sequences.fasta \
  >  --sequence-index filter/data/sequence_index.tsv \
  >  --metadata filter/data/metadata.tsv \
  >  --min-date 2012 \
  >  --group-by country year month \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 314159 \
  >  --output-strains "$TMP/filtered_strains_default.txt" > /dev/null
  WARNING: Asked to provide at most 5 sequences, but there are 8 groups.

By setting the subsample seed above, we should get the same results for both runs.

  $ diff -u <(sort "$TMP/filtered_strains_probabilistic.txt") <(sort "$TMP/filtered_strains_default.txt")
  $ rm -f "$TMP/filtered_strains_probabilistic.txt" "$TMP/filtered_strains_default.txt"
