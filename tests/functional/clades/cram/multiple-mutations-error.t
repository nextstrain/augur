Integration tests for augur clades.

  $ source "$TESTDIR"/_setup.sh

Multiple mutations at the same position on a single branch are a fatal error

  $ ${AUGUR} clades \
  >   --tree "$TESTDIR/../data/toy_tree.nwk" \
  >   --mutations "$TESTDIR/../data/toy_muts_multiple.json" \
  >   --clades "$TESTDIR/../data/toy_clades_nuc.tsv"
  ERROR: Multiple mutations at the same position on a single branch were found: Node A (nuc), Node AB (geneName)
  [2]
