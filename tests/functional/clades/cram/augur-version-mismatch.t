Integration tests for augur clades.

  $ source "$TESTDIR"/_setup.sh

Node-data JSONs produced from a different major version of augur
are not allowed.

  $ ${AUGUR} clades \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --mutations "$TESTDIR/../data/aa_muts_generated_by.json" \
  >   --clades "$TESTDIR/../data/clades.tsv" \
  >   --output-node-data clades.json
  ERROR: Augur version incompatibility detected: the JSON .*aa_muts_generated_by\.json.* was generated by \{'program': 'augur', 'version': '21.1.0'\}, which is incompatible with the current augur version \([.0-9]+\). We suggest you rerun the pipeline using the current version of augur. (re)
  [2]

Skipping validation allows mismatched augur versions to be used without error.

  $ ${AUGUR} clades \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --mutations "$TESTDIR/../data/aa_muts_generated_by.json" \
  >   --clades "$TESTDIR/../data/clades.tsv" \
  >   --output-node-data clades.json \
  >   --skip-validation &>/dev/null
