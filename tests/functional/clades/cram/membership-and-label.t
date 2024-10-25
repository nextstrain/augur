Integration tests for augur clades.

  $ source "$TESTDIR"/_setup.sh

Test custom membership key + label key. The only change should be the key names

  $ ${AUGUR} clades \
  >   --tree "$TESTDIR/../data/tree.nwk" \
  >   --mutations "$TESTDIR/../data/aa_muts.json" "$TESTDIR/../data/nt_muts_small.json" \
  >   --clades "$TESTDIR/../data/clades.tsv" \
  >   --membership-name lineage --label-name origin \
  >   --output-node-data clades_custom.json &>/dev/null

  $ cat clades_custom.json | sed "s/lineage/clade_membership/" | sed "s/origin/clade/" > clades_sed.json

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/clades.json" clades_sed.json \
  >   --exclude-paths "root['generated_by']"
  {}
