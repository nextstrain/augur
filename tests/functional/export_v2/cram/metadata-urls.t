Setup

  $ source "$TESTDIR"/_setup.sh

Create files for testing.

  $ cat >metadata.tsv <<~~
  > strain	field_A
  > tipA	nextstrain: https://nextstrain.org
  > tipB	BB: not-a-url
  > tipC	github <https://github.com>
  > tipD	DD <not-a-url>
  > tipE	EE
  > tipF	FF
  > ~~

  $ cat >tree.nwk <<~~
  > (tipA:1,(tipB:1,tipC:1)internalBC:2,(tipD:3,tipE:4,tipF:1)internalDEF:5)ROOT:0;
  > ~~

Check that URLs were extracted from metadata values when added as an "extra metadata column"

  $ ${AUGUR} export v2 \
  >  --tree tree.nwk \
  >  --metadata metadata.tsv \
  >  --metadata-columns "field_A" \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset.json
  Validating produced JSON
  Validating schema of 'dataset.json'...
  Validating that the JSON is internally consistent...
  Validation of 'dataset.json' succeeded.
  

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-with-parsed-urls.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']"
  {}


Check that URLs were extracted from metadata values when used as a coloring

  $ ${AUGUR} export v2 \
  >  --tree tree.nwk \
  >  --metadata metadata.tsv \
  >  --color-by-metadata "field_A" \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset2.json
  Trait 'field_A' was guessed as being type 'categorical'. Use a 'config' file if you'd like to set this yourself.
  Validating produced JSON
  Validating schema of 'dataset2.json'...
  Validating that the JSON is internally consistent...
  Validation of 'dataset2.json' succeeded.
  

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-with-parsed-urls.json" dataset2.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['colorings']" "root['meta']['filters']"
  {}