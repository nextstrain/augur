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
    .meta.updated "invalid-date" failed pattern validation for "^[0-9X]{4}-[0-9X]{2}-[0-9X]{2}$"
  FATAL ERROR: Validation of 'bad-date.json' failed.
  [2]
