Integration tests for augur ancestral.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="${AUGUR:-../../bin/augur}"

Infer ancestral sequences for the given tree and alignment.
The default is to infer ambiguous bases, so there should not be N bases in the inferred output sequences.

  $ ${AUGUR} ancestral \
  >  --tree ancestral/tree_raw.nwk \
  >  --alignment ancestral/aligned.fasta \
  >  --output-node-data "$TMP/ancestral_mutations.json" \
  >  --output-sequences "$TMP/ancestral_sequences.fasta" > /dev/null

  $ grep N "$TMP/ancestral_sequences.fasta"
  >NODE_0000000

Infer ancestral sequences for the given tree and alignment, explicitly requesting that ambiguous bases are inferred.
There should not be N bases in the inferred output sequences.

  $ ${AUGUR} ancestral \
  >  --tree ancestral/tree_raw.nwk \
  >  --alignment ancestral/aligned.fasta \
  >  --infer-ambiguous \
  >  --output-node-data "$TMP/ancestral_mutations.json" \
  >  --output-sequences "$TMP/ancestral_sequences.fasta" > /dev/null

  $ grep N "$TMP/ancestral_sequences.fasta"
  >NODE_0000000

Infer ancestral sequences for the given tree and alignment, explicitly requesting that ambiguous bases are NOT inferred.
There be N bases in the inferred output sequences.

  $ ${AUGUR} ancestral \
  >  --tree ancestral/tree_raw.nwk \
  >  --alignment ancestral/aligned.fasta \
  >  --keep-ambiguous \
  >  --output-node-data "$TMP/ancestral_mutations.json" \
  >  --output-sequences "$TMP/ancestral_sequences.fasta" > /dev/null

  $ grep N "$TMP/ancestral_sequences.fasta"
  >NODE_0000000
  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGTATCAACAGGTTTT
  AGCTGTGGATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
  NNNNNNNNNNNNNNNNNNNNNNNNNNNNN
  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGTATCAACAGGTTTT
  AGCTGTGGANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
  NNNNNNNNNNNNNNNNNNNNNNNNNNNNN

  $ popd > /dev/null
