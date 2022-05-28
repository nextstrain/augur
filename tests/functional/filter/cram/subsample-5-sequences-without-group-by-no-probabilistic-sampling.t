Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Filter with subsampling where no more than 5 sequences are requested and no groups are specified.
This generates a dummy category and subsamples from there. With no-probabilistic-sampling we expect exactly 5 sequences.

Pandas engine
-------------

  $ ${AUGUR} filter \
  >  --sequences filter/data/sequences.fasta \
  >  --sequence-index filter/data/sequence_index.tsv \
  >  --metadata filter/data/metadata.tsv \
  >  --min-date 2012 \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 314159 \
  >  --no-probabilistic-sampling \
  >  --output "$TMP/filtered.fasta" > /dev/null
  $ grep ">" "$TMP/filtered.fasta" | wc -l
  \s*5 (re)
  $ rm -f "$TMP/filtered.fasta"

SQLite engine
-------------

  $ ${AUGUR} filter --engine sqlite \
  >  --sequences filter/data/sequences.fasta \
  >  --sequence-index filter/data/sequence_index.tsv \
  >  --metadata filter/data/metadata.tsv \
  >  --min-date 2012 \
  >  --subsample-max-sequences 5 \
  >  --subsample-seed 314159 \
  >  --no-probabilistic-sampling \
  >  --output "$TMP/filtered.fasta" > /dev/null
  $ grep ">" "$TMP/filtered.fasta" | wc -l
  \s*5 (re)
  $ rm -f "$TMP/filtered.fasta"
