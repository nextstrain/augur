Setup

  $ source "$TESTDIR"/_setup.sh

Infer ancestral sequences for the given tree and alignment, explicitly requesting that ambiguous bases are inferred.
There should not be N bases in the inferred output sequences.

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --alignment $TESTDIR/../data/aligned.fasta \
  >  --infer-ambiguous \
  >  --seed 314159 \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations.json" \
  >  --output-sequences "$CRAMTMP/$TESTFILE/ancestral_sequences.fasta" > /dev/null

  $ grep "^N" "$CRAMTMP/$TESTFILE/ancestral_sequences.fasta"
  [1]
