Integration tests for augur clades2.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="${AUGUR:-../../bin/augur}"

Infer clades

  $ ${AUGUR} clades2 \
  >  --metadata "clades2/metadata.tsv" \
  >  --tree "clades2/tree.nwk" \
  >  --clade-column clade \
  >  --output-node-data "$TMP/clades.json" > /dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" "clades2/clades.json" "$TMP/clades.json"
  {}
  $ rm -f "$TMP/clades.json"

Infer clades where some tips are not in metadata

  $ ${AUGUR} clades2 \
  >  --metadata "clades2/partial_metadata.tsv" \
  >  --tree "clades2/tree.nwk" \
  >  --clade-column clade \
  >  --output-node-data "$TMP/clades.json" > /dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" "clades2/clades.json" "$TMP/clades.json"
  {}
  $ rm -f "$TMP/clades.json"
