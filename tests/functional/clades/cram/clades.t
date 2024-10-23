Integration tests for augur clades.

  $ source "$TESTDIR"/_setup.sh
  $ pushd "$TESTDIR" > /dev/null

Test augur clades with simple Zika input files and hierarchical clades.

  $ ${AUGUR} clades \
  >   --tree ../data/tree.nwk \
  >   --mutations ../data/aa_muts.json ../data/nt_muts_small.json \
  >   --clades ../data/clades.tsv \
  >   --output-node-data "$TMP/clades.json" &>/dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py"  ../data/clades.json "$TMP/clades.json" \
  >   --exclude-paths "root['generated_by']"
  {}

Test custom membership key + label key. The only change should be the key names

  $ ${AUGUR} clades \
  >   --tree ../data/tree.nwk \
  >   --mutations ../data/aa_muts.json ../data/nt_muts_small.json \
  >   --clades ../data/clades.tsv \
  >   --membership-name lineage --label-name origin \
  >   --output-node-data "$TMP/clades_custom.json" &>/dev/null

  $ cat "$TMP/clades_custom.json" | sed "s/lineage/clade_membership/" | sed "s/origin/clade/" > "$TMP/clades_sed.json"

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" ../data/clades.json "$TMP/clades_sed.json" \
  >   --exclude-paths "root['generated_by']"
  {}

Test the ability to _not_ export a branch label (same logic as not exporting the membership)

  $ ${AUGUR} clades \
  >   --tree ../data/tree.nwk \
  >   --mutations ../data/aa_muts.json ../data/nt_muts_small.json \
  >   --clades ../data/clades.tsv \
  >   --label-name none \
  >   --output-node-data "$TMP/clades_no-labels.json" &>/dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" ../data/clades.json "$TMP/clades_no-labels.json" \
  >   --exclude-paths "root['generated_by']"
  {'dictionary_item_removed': [root['branches']]}


A clade which exists at the root is not identified by inferring the root state
(i.e. we don't infer the root state to be A if we observe a subsequent A10T mutation)
This is an oversight and ideally would be fixed

  $ ${AUGUR} clades \
  >   --tree ../data/toy_tree.nwk \
  >   --mutations ../data/toy_muts_no_ref.json \
  >   --clades ../data/toy_clades_nuc.tsv \
  >   --output-node-data "$TMP/toy_clades_1.json" &>/dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" ../data/toy_clades_1.json "$TMP/toy_clades_1.json" \
  >   --exclude-paths "root['generated_by']"
  {}

A clade which exists at the root is identified (and correctly propogated) if the root sequence
is explicitly set.

  $ ${AUGUR} clades \
  >   --tree ../data/toy_tree.nwk \
  >   --mutations ../data/toy_muts_ref.json \
  >   --clades ../data/toy_clades_nuc.tsv \
  >   --output-node-data "$TMP/toy_clades_2a.json" &>/dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" ../data/toy_clades_2.json "$TMP/toy_clades_2a.json" \
  >   --exclude-paths "root['generated_by']"
  {}

A clade which exists at the root is identified (and correctly propogated) without a root sequence
if the (branch leading to the) root has the clade-defining mutation.

  $ ${AUGUR} clades \
  >   --tree ../data/toy_tree.nwk \
  >   --mutations ../data/toy_muts_explicit_root_mutation.json \
  >   --clades ../data/toy_clades_nuc.tsv \
  >   --output-node-data "$TMP/toy_clades_2b.json" &>/dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" ../data/toy_clades_2.json "$TMP/toy_clades_2b.json" \
  >   --exclude-paths "root['generated_by']"
  {}

Multiple mutations at the same position on a single branch are a fatal error

  $ ${AUGUR} clades \
  >   --tree ../data/toy_tree.nwk \
  >   --mutations ../data/toy_muts_multiple.json \
  >   --clades ../data/toy_clades_nuc.tsv
  ERROR: Multiple mutations at the same position on a single branch were found: Node A (nuc), Node AB (geneName)
  [2]
