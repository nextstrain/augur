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