Integration tests for augur tree.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="../../bin/augur"

Try building a tree with IQ-TREE.

  $ ${AUGUR} tree \
  >  --alignment tree/aligned.fasta \
  >  --method iqtree \
  >  --output "$TMP/tree_raw.nwk" \
  >  --nthreads 1 > /dev/null

Try building a tree with IQ-TREE with more threads (4) than there are input sequences (3).

  $ ${AUGUR} tree \
  >  --alignment tree/aligned.fasta \
  >  --method iqtree \
  >  --output "$TMP/tree_raw.nwk" \
  >  --nthreads 4 > /dev/null
  WARNING: more threads requested than there are sequences; falling back to IQ-TREE's `-nt AUTO` mode.

  $ popd > /dev/null
