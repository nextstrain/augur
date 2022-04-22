Integration tests for augur index.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="../../bin/augur"

Index Zika sequences.

  $ ${AUGUR} index \
  >   --sequences index/sequences.fasta \
  >   --output "$TMP/sequence_index.tsv"

  $ diff -u "index/sequence_index.tsv" "$TMP/sequence_index.tsv"
  $ rm -f "$TMP/sequence_index.tsv"

Try indexing sequences that do not exist.
This should fail.

  $ ${AUGUR} index \
  >   --sequences index/missing_sequences.fasta \
  >   --output "$TMP/sequence_index.tsv"
  ERROR: Could not open sequences file 'index/missing_sequences.fasta'.
  [1]

  $ popd > /dev/null
