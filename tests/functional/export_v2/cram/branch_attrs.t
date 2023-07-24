Setup

  $ source "$TESTDIR"/_setup.sh

Test that attributes are correctly exported as branch_attrs. Currently this includes branch labels (node_data→branches),
mutations (node_data→nodes) and a historical node_data→nodes→<name>→clade_annotation branch label.

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --node-data "$TESTDIR/../data/div_node-data.json" "$TESTDIR/../data/nt_muts_1.json" "$TESTDIR/../data/aa_muts_1.json" "$TESTDIR/../data/branch-labels.json" \
  >   --maintainers "Nextstrain Team" \
  >   --output dataset.json > /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py"  "$TESTDIR/../data/dataset-with-branch-labels.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']"
  {}
