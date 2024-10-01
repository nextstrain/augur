Integration tests for augur index.

  $ source "$TESTDIR"/_setup.sh

Index Zika sequences.

  $ ${AUGUR} index \
  >   --sequences "$TESTDIR/index/sequences.fasta" \
  >   --output sequence_index.tsv

  $ diff -u "$TESTDIR/index/sequence_index.tsv" sequence_index.tsv
  $ rm -f sequence_index.tsv

Try indexing sequences that do not exist.
This should fail.

  $ ${AUGUR} index \
  >   --sequences index/missing_sequences.fasta \
  >   --output sequence_index.tsv
  ERROR: No such file or directory: 'index/missing_sequences.fasta'
  [2]

Try writing output to a directory that does not exist.
This should fail.

  $ ${AUGUR} index \
  >   --sequences "$TESTDIR/index/sequences.fasta" \
  >   --output results/sequence_index.tsv
  ERROR: No such file or directory: 'results/sequence_index.tsv'
  [2]
