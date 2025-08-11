Setup

  $ source "$TESTDIR"/_setup.sh

Create files for testing.

  $ cat >metadata.tsv <<~~
  > strain	field_A	field_A__url
  > tipA	nextstrain	https://nextstrain.org
  > tipB	BB: not-a-url	
  > tipC	github	https://github.com
  > tipD	DD <not-a-url>	
  > tipE	EE	invalid-url
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
  Validating schema of 'dataset.json'...
  Validation of 'dataset.json' succeeded.
  Validating produced JSON
  Validating that the JSON is internally consistent...
  

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
  Validating schema of 'dataset2.json'...
  Validation of 'dataset2.json' succeeded.
  Trait 'field_A' was guessed as being type 'categorical'. Use a 'config' file if you'd like to set this yourself.
  Validating produced JSON
  Validating that the JSON is internally consistent...
  

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-with-parsed-urls.json" dataset2.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['colorings']" "root['meta']['filters']"
  {}

Check that node-data JSONs, where non-special-case values are automatically used, are correctly handled.
The data is essentially the same, but tipB & tipD have empty-string URLs and tipF is missing the URL key entirely.

  $ cat >node-data.json <<~~
  > {"nodes":
  >   {
  >     "tipA": {"field_A": "nextstrain", "field_A__url": "https://nextstrain.org"},
  >     "tipB": {"field_A": "BB: not-a-url", "field_A__url": ""},
  >     "tipC": {"field_A": "github", "field_A__url": "https://github.com"},
  >     "tipD": {"field_A": "DD <not-a-url>", "field_A__url": ""},
  >     "tipE": {"field_A": "EE", "field_A__url": "invalid-url"},
  >     "tipF": {"field_A": "FF"}
  >   }
  > }
  > ~~

  $ ${AUGUR} export v2 \
  >  --tree tree.nwk \
  >  --node-data node-data.json \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset3.json
  Validating schema of 'dataset3.json'...
  Validation of 'dataset3.json' succeeded.
  Trait 'field_A' was guessed as being type 'categorical'. Use a 'config' file if you'd like to set this yourself.
  Validating produced JSON
  Validating that the JSON is internally consistent...
  


  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-with-parsed-urls.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']"
  {}