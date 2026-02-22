Setup

  $ source "$TESTDIR"/_setup.sh

Test a dataset with only AA-mutations and no `nuc` annotations block, i.e. a protein analysis

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --node-data "$TESTDIR/../data/div_node-data.json" "$TESTDIR/../data/aa_muts_2.json" \
  >   --maintainers "Nextstrain Team" \
  >   --output dataset.json &>/dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py"  "$TESTDIR/../data/dataset-with-only-amino-acids.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']"
  {}

If we don't include the `nuc` block but do include nuc mutations we should show a warning

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --node-data "$TESTDIR/../data/div_node-data.json" "$TESTDIR/../data/nt_muts_2.json" \
  >   --maintainers "Nextstrain Team" \
  >   --output dataset2.json 1>/dev/null 2>dataset2.log

  $ grep -q "WARNING:  The tree defined mutations on gene nuc which doesn't appear in the metadata annotations object." dataset2.log || \
  >  echo "Test failure: export v2 should have warned about a missing 'nuc' annotation block"
