Setup

  $ source "$TESTDIR"/_setup.sh

Test date pattern validation.

  $ cat >bad-date.json <<~~
  > {
  >   "version": "v2",
  >   "meta": {
  >     "updated": "invalid-date",
  >     "panels": ["tree"]
  >   },
  >   "tree": {
  >     "name": "test_tree",
  >     "node_attrs": {"div": 0}
  >   }
  > }
  > ~~

  $ ${AUGUR} validate export-v2 bad-date.json
  Validating schema of 'bad-date.json'...
    .meta.updated failed: Expected a date in the format YYYY-MM-DD but found "invalid-date"
  FATAL ERROR: Validation of 'bad-date.json' failed.
  [2]
