Setup

  $ source "$TESTDIR"/_setup.sh

Test malformed JSON input:

  $ cat >malformed.json <<~~
  > {
  >   "version": "v2",
  >   "meta": {
  >     "updated": "2023-12-25",
  >     "panels": ["tree"]
  >   },
  >   "tree": {
  >     "name": "test_tree",
  >     "node_attrs": {"div": 0}
  >   }
  > ~~

  $ ${AUGUR} validate export-v2 malformed.json
  FATAL ERROR: Supplied JSON to validate (malformed.json) is not a valid JSON
  [2]

Test JSON with trailing comma:

  $ cat >trailing-comma.json <<~~
  > {
  >   "version": "v2",
  >   "meta": {
  >     "updated": "2023-12-25",
  >     "panels": ["tree"],
  >   },
  >   "tree": {
  >     "name": "test_tree",
  >     "node_attrs": {"div": 0}
  >   }
  > }
  > ~~

  $ ${AUGUR} validate export-v2 trailing-comma.json
  FATAL ERROR: Supplied JSON to validate (trailing-comma.json) is not a valid JSON
  [2]

Test completely invalid JSON:

  $ cat >invalid.json <<~~
  > not json at all
  > ~~

  $ ${AUGUR} validate export-v2 invalid.json
  FATAL ERROR: Supplied JSON to validate (invalid.json) is not a valid JSON
  [2]

Test JSON with duplicate keys:

  $ cat >duplicate-keys.json <<~~
  > {
  >   "version": "v2",
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

  $ ${AUGUR} validate export-v2 duplicate-keys.json
  Validating schema of 'duplicate-keys.json'...
    .version failed: Expected 'v2' but found 'v1'
  FATAL ERROR: Validation of 'duplicate-keys.json' failed.
  [2]
