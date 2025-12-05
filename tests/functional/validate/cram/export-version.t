Setup

  $ source "$TESTDIR"/_setup.sh

Test with incorrect version.

  $ cat >wrong-version.json <<~~
  > {
  >   "version": "v1",
  >   "meta": {
  >     "updated": "2023-12-25",
  >     "panels": ["tree"]
  >   },
  >   "tree": {
  >     "name": "test_tree",
  >     "node_attrs": {"div": 0}
  >   }
  > }
  > ~~

  $ ${AUGUR} validate export-v2 wrong-version.json
  Validating schema of 'wrong-version.json'...
    .version failed: Expected 'v2' but found 'v1'
  FATAL ERROR: Validation of 'wrong-version.json' failed.
  [2]
