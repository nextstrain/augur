Integration tests for augur traits.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="../../bin/augur"

Infer the ancestral region for a given tree and metadata.

  $ ${AUGUR} traits \
  >  --metadata "traits/metadata.tsv" \
  >  --tree "traits/tree.nwk" \
  >  --columns region \
  >  --output-node-data "$TMP/traits.json" > /dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" "traits/traits_region.json" "$TMP/traits.json" --significant-digits 5
  {}
  $ rm -f "$TMP/traits.json"

Infer the ancestral "virus" value from the same metadata.
Since there is only a single virus in the data, Augur warns the user through stderr.

  $ ${AUGUR} traits \
  >  --metadata "traits/metadata.tsv" \
  >  --tree "traits/tree.nwk" \
  >  --columns virus \
  >  --output-node-data "$TMP/traits.json" > /dev/null
  WARNING: only one state found for discrete state reconstruction: ['zika']

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" "traits/traits_virus.json" "$TMP/traits.json" --significant-digits 5
  {}
  $ rm -f "$TMP/traits.json"

Repeat inference of a trait with a single value, but request confidence intervals.
This should similarly warn the user through stderr, but it should produce an error.

  $ ${AUGUR} traits \
  >  --metadata "traits/metadata.tsv" \
  >  --tree "traits/tree.nwk" \
  >  --columns virus \
  >  --confidence \
  >  --output-node-data "$TMP/traits.json" > /dev/null
  WARNING: only one state found for discrete state reconstruction: ['zika']

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" "traits/traits_virus.json" "$TMP/traits.json" --significant-digits 5
  {}
  $ rm -f "$TMP/traits.json"

Switch back to the original directory where testing started.

  $ popd > /dev/null
