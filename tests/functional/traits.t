Integration tests for augur traits.

  $ source "$TESTDIR"/_setup.sh

Infer the ancestral region for a given tree and metadata.

  $ ${AUGUR} traits \
  >  --metadata "$TESTDIR/traits/metadata.tsv" \
  >  --tree "$TESTDIR/traits/tree.nwk" \
  >  --columns region \
  >  --output-node-data "traits.json" > /dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" "$TESTDIR/traits/traits_region.json" "traits.json" --significant-digits 5
  {}

Infer the ancestral region for a tree and metadata where one or more records have a missing value ("?") in the region field.
Tips with missing values should get their values inferred, too.
In this case, a sample from Panama (North America) has its region inferred as "South America".

  $ ${AUGUR} traits \
  >  --metadata "$TESTDIR/traits/metadata_with_missing_region.tsv" \
  >  --tree "$TESTDIR/traits/tree.nwk" \
  >  --columns region \
  >  --output-node-data "traits.json" > /dev/null

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" "$TESTDIR/traits/traits_with_missing_region.json" "traits.json" --significant-digits 2
  {}

Infer the ancestral "virus" value from the same metadata.
Since there is only a single virus in the data, Augur warns the user through stderr.

  $ ${AUGUR} traits \
  >  --metadata "$TESTDIR/traits/metadata.tsv" \
  >  --tree "$TESTDIR/traits/tree.nwk" \
  >  --columns virus \
  >  --output-node-data "traits.json" > /dev/null
  WARNING: only one state found for discrete state reconstruction: ['zika']

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" "$TESTDIR/traits/traits_virus.json" "traits.json" --significant-digits 5
  {}

Repeat inference of a trait with a single value, but request confidence intervals.
This should similarly warn the user through stderr, but it should produce an error.

  $ ${AUGUR} traits \
  >  --metadata "$TESTDIR/traits/metadata.tsv" \
  >  --tree "$TESTDIR/traits/tree.nwk" \
  >  --columns virus \
  >  --confidence \
  >  --output-node-data "traits.json" > /dev/null
  WARNING: only one state found for discrete state reconstruction: ['zika']

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" "$TESTDIR/traits/traits_virus.json" "traits.json" --significant-digits 5
  {}

Infer the ancestral "virus" value from the metadata after replacing the "zika" values with missing data values ("?").
Augur should warn that there were no discrete states found for reconstruction, since "?" is not a valid state on its own.

  $ sed 's/zika/?/' "$TESTDIR/traits/metadata.tsv" > "metadata_with_missing_virus.tsv"
  $ ${AUGUR} traits \
  >  --metadata "metadata_with_missing_virus.tsv" \
  >  --tree "$TESTDIR/traits/tree.nwk" \
  >  --columns virus \
  >  --output-node-data "traits.json" > /dev/null
  WARNING: no states found for discrete state reconstruction.
