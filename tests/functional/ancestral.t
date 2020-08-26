Integration tests for augur ancestral.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="../../bin/augur"

Infer ancestral sequences for the given tree and alignment.
The default is to infer ambiguous bases, so there should not be N bases in the inferred output sequences.

  $ ${AUGUR} ancestral \
  >  --tree ancestral/tree_raw.nwk \
  >  --alignment ancestral/aligned.fasta \
  >  --output-node-data "$TMP/ancestral_mutations.json" \
  >  --output-sequences "$TMP/ancestral_sequences.fasta" > /dev/null
  */treetime/aa_models.py:108: VisibleDeprecationWarning: Creating an ndarray from ragged nested sequences (which is a list-or-tuple of lists-or-tuples-or ndarrays with different lengths or shapes) is deprecated. If you meant to do this, you must specify 'dtype=object' when creating the ndarray (glob)
    [0.800038509951648, 1.20274751778601, 1.55207513886163, 1.46600946033173, 0.830022143283238, 1.5416250309563, 1.53255698189437, 1.41208067821187, 1.47469999960758, 0.351200119909572, 0.570542199221932, 1.21378822764856, 0.609532859331199, 0.692733248746636, 1.40887880416009, 1.02015839286433, 0.807404666228614, 1.268589159299, 0.933095433689795]

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
  */treetime/aa_models.py:108: VisibleDeprecationWarning: Creating an ndarray from ragged nested sequences (which is a list-or-tuple of lists-or-tuples-or ndarrays with different lengths or shapes) is deprecated. If you meant to do this, you must specify 'dtype=object' when creating the ndarray (glob)
    [0.800038509951648, 1.20274751778601, 1.55207513886163, 1.46600946033173, 0.830022143283238, 1.5416250309563, 1.53255698189437, 1.41208067821187, 1.47469999960758, 0.351200119909572, 0.570542199221932, 1.21378822764856, 0.609532859331199, 0.692733248746636, 1.40887880416009, 1.02015839286433, 0.807404666228614, 1.268589159299, 0.933095433689795]

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
  */treetime/aa_models.py:108: VisibleDeprecationWarning: Creating an ndarray from ragged nested sequences (which is a list-or-tuple of lists-or-tuples-or ndarrays with different lengths or shapes) is deprecated. If you meant to do this, you must specify 'dtype=object' when creating the ndarray (glob)
    [0.800038509951648, 1.20274751778601, 1.55207513886163, 1.46600946033173, 0.830022143283238, 1.5416250309563, 1.53255698189437, 1.41208067821187, 1.47469999960758, 0.351200119909572, 0.570542199221932, 1.21378822764856, 0.609532859331199, 0.692733248746636, 1.40887880416009, 1.02015839286433, 0.807404666228614, 1.268589159299, 0.933095433689795]

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
