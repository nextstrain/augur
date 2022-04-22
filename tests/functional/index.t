Integration tests for augur index.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="../../bin/augur"

Index Zika sequences.

  $ ${AUGUR} index \
  >   --sequences index/sequences.fasta \
  >   --output "$TMP/sequence_index.tsv"

  $ diff -u -w "index/sequence_index.tsv" "$TMP/sequence_index.tsv"
  $ rm -f "$TMP/sequence_index.tsv"

Index sequences with duplicate strains.
We should only include the first instance of any duplicate in the output and should warn the user about duplicates.

  $ ${AUGUR} index \
  >   --sequences index/sequences_with_duplicate.fasta \
  >   --output "$TMP/sequence_index.tsv"
  WARNING: Found duplicate sequence for strain 'PRVABC59'. Only the first sequence for this strain will be indexed.

  $ diff -u -w "index/sequence_index.tsv" "$TMP/sequence_index.tsv"
  $ rm -f "$TMP/sequence_index.tsv"

  $ popd > /dev/null
