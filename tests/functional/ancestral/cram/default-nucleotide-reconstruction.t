Setup

  $ source "$TESTDIR"/_setup.sh

Infer ancestral sequences for the given tree and alignment.
The default is to infer ambiguous bases, so there should not be N bases in the inferred output sequences.

  $ ${AUGUR} ancestral \
  >  --tree $TESTDIR/../data/tree.nwk \
  >  --alignment $TESTDIR/../data/aligned.fasta \
  >  --seed 314159 \
  >  --output-node-data "$CRAMTMP/$TESTFILE/ancestral_mutations.json" \
  >  --output-sequences "$CRAMTMP/$TESTFILE/ancestral_sequences.fasta" > /dev/null

  $ grep "^N" "$CRAMTMP/$TESTFILE/ancestral_sequences.fasta"
  [1]

Check that the reference length was correctly exported as the nuc annotation

  $ grep -A 6 'annotations' "$CRAMTMP/$TESTFILE/ancestral_mutations.json"
    "annotations": {
      "nuc": {
        "end": 10769,
        "start": 1,
        "strand": "+",
        "type": "source"
      }
