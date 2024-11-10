Integration tests for augur clades.

  $ source "$TESTDIR"/_setup.sh

Test the ability to _not_ export a branch label (same logic as not exporting the membership)

  $ ${AUGUR} clades \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --mutations "$TESTDIR/../data/aa_muts.json" "$TESTDIR/../data/nt_muts_small.json" \
  >   --clades "$TESTDIR/../data/clades.tsv" \
  >   --label-name none \
  >   --output-node-data clades_no-labels.json &>/dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/clades.json" clades_no-labels.json \
  >   --exclude-paths "root['generated_by']"
  {'dictionary_item_removed': [root['branches']]}
