Setup

  $ source "$TESTDIR"/_setup.sh

Create files for testing.

  $ cat >metadata.tsv <<~~
  > strain	field_A	field_A__url	field_B	field_B__url
  > tipA	nextstrain	https://nextstrain.org	1	https://nextstrain.org
  > tipB	BB: not-a-url		2	
  > tipC	github	https://github.com	3	https://github.com
  > tipD	DD <not-a-url>		4	
  > tipE	EE	invalid-url	5	invalid-url
  > tipF	FF		6	
  > ~~

  $ cat >tree.nwk <<~~
  > (tipA:1,(tipB:1,tipC:1)internalBC:2,(tipD:3,tipE:4,tipF:1)internalDEF:5)ROOT:0;
  > ~~

Check that URLs were extracted from metadata values when added as an "extra metadata column"

  $ ${AUGUR} export v2 \
  >  --tree tree.nwk \
  >  --metadata metadata.tsv \
  >  --metadata-columns "field_A" "field_B" \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset.json &> /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-with-parsed-urls.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']"
  {}


Check that URLs were extracted from metadata values when used as a coloring

  $ ${AUGUR} export v2 \
  >  --tree tree.nwk \
  >  --metadata metadata.tsv \
  >  --color-by-metadata "field_A" "field_B" \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset2.json &> /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-with-parsed-urls.json" dataset2.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['colorings']" "root['meta']['filters']"
  {}

Check that node-data JSONs, where non-special-case values are automatically used, are correctly handled.
The data is essentially the same, but tipB & tipD have empty-string URLs and tipF is missing the URL key entirely.

  $ cat >node-data.json <<~~
  > {"nodes":
  >   {
  >     "tipA": {"field_A": "nextstrain", "field_A__url": "https://nextstrain.org", "field_B": 1, "field_B__url": "https://nextstrain.org"},
  >     "tipB": {"field_A": "BB: not-a-url", "field_A__url": "", "field_B": 2, "field_B__url": ""},
  >     "tipC": {"field_A": "github", "field_A__url": "https://github.com", "field_B": 3, "field_B__url": "https://github.com"},
  >     "tipD": {"field_A": "DD <not-a-url>", "field_A__url": "", "field_B": 4, "field_B__url": ""},
  >     "tipE": {"field_A": "EE", "field_A__url": "invalid-url", "field_B": 5, "field_B__url": "invalid-url"},
  >     "tipF": {"field_A": "FF", "field_B": 6}
  >   }
  > }
  > ~~

  $ ${AUGUR} export v2 \
  >  --tree tree.nwk \
  >  --node-data node-data.json \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset3.json &> /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-with-parsed-urls.json" dataset3.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['colorings']" "root['meta']['filters']"
  {}
