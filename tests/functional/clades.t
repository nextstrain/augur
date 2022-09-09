Integration tests for augur clades.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="${AUGUR:-../../bin/augur}"

Test augur clades with simple Zika input files and hierarchical clades.

  $ ${AUGUR} clades \
  >   --tree clades/tree.nwk \
  >   --mutations clades/aa_muts.json clades/nt_muts_small.json \
  >   --clades clades/clades.tsv \
  >   --output-node-data "$TMP/clades.json" &>/dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py"  clades/clades.json "$TMP/clades.json" \
  >   --exclude-paths "root['generated_by']"
  {}

Test custom membership key + label key. The only change should be the key names

  $ ${AUGUR} clades \
  >   --tree clades/tree.nwk \
  >   --mutations clades/aa_muts.json clades/nt_muts_small.json \
  >   --clades clades/clades.tsv \
  >   --membership-name lineage --label-name origin \
  >   --output-node-data "$TMP/clades_custom.json" &>/dev/null

  $ cat "$TMP/clades_custom.json" | sed "s/lineage/clade_membership/" | sed "s/origin/clade/" > "$TMP/clades_sed.json"

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" clades/clades.json "$TMP/clades_sed.json" \
  >   --exclude-paths "root['generated_by']"
  {}

Test the ability to _not_ export a branch label (same logic as not exporting the membership)

  $ ${AUGUR} clades \
  >   --tree clades/tree.nwk \
  >   --mutations clades/aa_muts.json clades/nt_muts_small.json \
  >   --clades clades/clades.tsv \
  >   --label-name none \
  >   --output-node-data "$TMP/clades_no-labels.json" &>/dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" clades/clades.json "$TMP/clades_no-labels.json" \
  >   --exclude-paths "root['generated_by']"
  {'dictionary_item_removed': [root['branches']]}
