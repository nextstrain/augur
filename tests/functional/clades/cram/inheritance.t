Integration tests for augur clades.

  $ source "$TESTDIR"/_setup.sh

  $ ${AUGUR} clades \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --mutations "$TESTDIR/../data/inheritance-muts.json" \
  >   --clades "$TESTDIR/../data/inheritance-clades.tsv" \
  >   --output-node-data clades.json &>/dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py"  "$TESTDIR/../data/inheritance-clades.json" clades.json \
  >   --exclude-paths "root['generated_by']"
  {}
