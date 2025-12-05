Setup

  $ source "$TESTDIR"/_setup.sh

Test with a missing closing }.

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

Test with a trailing comma.

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

Test completely invalid JSON.

  $ cat >invalid.json <<~~
  > not json at all
  > ~~

  $ ${AUGUR} validate export-v2 invalid.json
  FATAL ERROR: Supplied JSON to validate (invalid.json) is not a valid JSON
  [2]

If there are two entries of the same name, the last one will be used.

  $ cat >duplicate-keys.json <<~~
  > {
  >   "version": "v1",
  >   "version": "v2",
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
  	WARNING:  No colourings were provided
  Validation of 'duplicate-keys.json' succeeded, but there were warnings you may want to resolve.
  Validating that the JSON is internally consistent...
