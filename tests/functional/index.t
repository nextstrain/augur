Integration tests for augur index.

  $ source "$TESTDIR"/_setup.sh

Index Zika sequences.

  $ ${AUGUR} index \
  >   --sequences "$TESTDIR/index/sequences.fasta" \
  >   --output sequence_index.tsv

  $ diff -u "$TESTDIR/index/sequence_index.tsv" sequence_index.tsv
  $ rm -f sequence_index.tsv

Index protein (HA1) sequences

  $ ${AUGUR} index \
  >   --seq-type aa \
  >   --sequences "$TESTDIR/index/HA1.fasta" \
  >   --output sequence_index_HA1.tsv

  $ diff -u "$TESTDIR/index/HA1_index.tsv" sequence_index_HA1.tsv
  $ rm -f sequence_index_HA1.tsv

Specifying a seq-type not 'aa' or 'nuc' fails
(The exact error message isn't written out below as the quoting was different on my system vs. CI)

  $ ${AUGUR} index \
  >   --seq-type protein \
  >   --sequences "$TESTDIR/index/HA1.fasta" \
  >   --output sequence_index_HA1.tsv
  usage: augur index .+ (re)
  .+ (re)
  augur index: error: argument --seq-type: invalid choice: 'protein' .+ (re)
  [2]

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


Empty sequences files are invalid

  $ touch empty.fasta

  $ ${AUGUR} index \
  >   --sequences empty.fasta \
  >   --output sequence_index_empty.tsv
  ERROR: Empty sequences provided
  [2]
