Integration tests for augur clades.

  $ source "$TESTDIR"/_setup.sh

A clade which exists at the root is not identified by inferring the root state
(i.e. we don't infer the root state to be A if we observe a subsequent A10T mutation)
This is an oversight and ideally would be fixed

  $ ${AUGUR} clades \
  >   --tree "$TESTDIR/../data/toy_tree.nwk" \
  >   --mutations "$TESTDIR/../data/toy_muts_no_ref.json" \
  >   --clades "$TESTDIR/../data/toy_clades_nuc.tsv" \
  >   --output-node-data toy_clades_1.json &>/dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/toy_clades_1.json" toy_clades_1.json \
  >   --exclude-paths "root['generated_by']"
  {}

A clade which exists at the root is identified (and correctly propogated) if the root sequence
is explicitly set.

  $ ${AUGUR} clades \
  >   --tree "$TESTDIR/../data/toy_tree.nwk" \
  >   --mutations "$TESTDIR/../data/toy_muts_ref.json" \
  >   --clades "$TESTDIR/../data/toy_clades_nuc.tsv" \
  >   --output-node-data toy_clades_2a.json &>/dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/toy_clades_2.json" toy_clades_2a.json \
  >   --exclude-paths "root['generated_by']"
  {}

A clade which exists at the root is identified (and correctly propogated) without a root sequence
if the (branch leading to the) root has the clade-defining mutation.

  $ ${AUGUR} clades \
  >   --tree "$TESTDIR/../data/toy_tree.nwk" \
  >   --mutations "$TESTDIR/../data/toy_muts_explicit_root_mutation.json" \
  >   --clades "$TESTDIR/../data/toy_clades_nuc.tsv" \
  >   --output-node-data toy_clades_2b.json &>/dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/toy_clades_2.json" toy_clades_2b.json \
  >   --exclude-paths "root['generated_by']"
  {}
