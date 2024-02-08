Setup

  $ source "$TESTDIR"/_setup.sh

Create files for testing.

  $ cat >metadata.tsv <<~~
  > strain	field_A	field_B
  > tipA	AA	AAA
  > tipB	BB	BBB
  > tipC	CC	CCC
  > tipD	DD	DDD
  > tipE	EE	EEE
  > tipF	FF	FFF
  > ~~

  $ cat >tree.nwk <<~~
  > (tipA:1,(tipB:1,tipC:1)internalBC:2,(tipD:3,tipE:4,tipF:1)internalDEF:5)ROOT:0;
  > ~~

  $ cat >auspice-config.json <<~~
  > {"metadata_columns": ["field_A", "field_B"]}
  > ~~

  $ cat >auspice-config-overridden.json <<~~
  > {"metadata_columns": ["overridden_field"]}
  > ~~

Run export with tree and metadata with additional columns.

  $ ${AUGUR} export v2 \
  >  --tree tree.nwk \
  >  --metadata metadata.tsv \
  >  --metadata-columns "field_A" "field_B" \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset.json > /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-with-additional-metadata-columns.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['maintainers']"
  {}

Missing columns are skipped with a warning.

  $ ${AUGUR} export v2 \
  >  --tree tree.nwk \
  >  --metadata metadata.tsv \
  >  --metadata-columns "field_A" "field_B" "missing_field" \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset.json > /dev/null
  WARNING: Requested metadata column 'missing_field' does not exist and will not be exported
  \s{0} (re)

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-with-additional-metadata-columns.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['maintainers']"
  {}

Specifying a fields with both --metadata-columns and --colory-by-metadata should result in field used as a coloring and a filter.

  $ ${AUGUR} export v2 \
  >  --tree tree.nwk \
  >  --metadata metadata.tsv \
  >  --metadata-columns "field_A" "field_B" \
  >  --color-by-metadata "field_B" \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset.json > /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-with-additional-metadata-columns.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['maintainers']"
  {'iterable_item_added': {"root['meta']['colorings'][0]": {'key': 'field_B', 'title': 'field_B', 'type': 'categorical'}, "root['meta']['filters'][0]": 'field_B'}}

Missing columns are skipped with a warning when specified by both --metadata-columns and --color-by-metadata.

  $ ${AUGUR} export v2 \
  >  --tree tree.nwk \
  >  --metadata metadata.tsv \
  >  --metadata-columns "field_A" "field_B" "missing_field" \
  >  --color-by-metadata "missing_field" \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset.json > /dev/null
  WARNING: Requested metadata column 'missing_field' does not exist and will not be exported
  \s{0} (re)
  WARNING: Requested color-by field 'missing_field' does not exist and will not be used as a coloring or exported.
  \s{0} (re)

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-with-additional-metadata-columns.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['maintainers']"
  {}

Specifying additional metadata columns via the Auspice configuration file.

  $ ${AUGUR} export v2 \
  >  --tree tree.nwk \
  >  --metadata metadata.tsv \
  >  --auspice-config auspice-config.json \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset.json > /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-with-additional-metadata-columns.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['maintainers']"
  {}

Specifying additional metadata columns via command line overrides the Auspice configuration file.

  $ ${AUGUR} export v2 \
  >  --tree tree.nwk \
  >  --metadata metadata.tsv \
  >  --auspice-config auspice-config-overridden.json \
  >  --metadata-columns "field_A" "field_B" \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset.json > /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-with-additional-metadata-columns.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['maintainers']"
  {}
