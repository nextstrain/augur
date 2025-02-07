Setup

  $ source "$TESTDIR"/_setup.sh

Create files for testing.

  $ cat >metadata.tsv <<~~
  > strain	field_A	field_B
  > 1.00	AA	AAA
  > 2.00	BB	BBB
  > 3.00	CC	CCC
  > 4.00	DD	DDD
  > 5.00	EE	EEE
  > 6.00	FF	FFF
  > ~~

  $ cat >tree.nwk <<~~
  > (1.00:1,(2.00:1,3.00:1)internalBC:2,(4.00:3,5.00:4,6.00:1)internalDEF:5)ROOT:0;
  > ~~

Run export with tree and metadata with additional columns.
The metadata should match with the tree even though the names are floats
because we force the index column to be strings.

  $ ${AUGUR} export v2 \
  >  --tree tree.nwk \
  >  --metadata metadata.tsv \
  >  --metadata-columns "field_A" "field_B" \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset.json > /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-with-float-strains.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['maintainers']"
  {}
