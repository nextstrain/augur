Integration tests for augur clades.

  $ source "$TESTDIR"/_setup.sh

Test augur clades with simple Zika input files and hierarchical clades.

  $ ${AUGUR} clades \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --mutations "$TESTDIR/../data/aa_muts.json" "$TESTDIR/../data/nt_muts_small.json" \
  >   --clades "$TESTDIR/../data/clades.tsv" \
  >   --output-node-data clades.json &>/dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py"  "$TESTDIR/../data/clades.json" clades.json \
  >   --exclude-paths "root['generated_by']"
  {}
