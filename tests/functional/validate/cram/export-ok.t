Setup

  $ source "$TESTDIR"/_setup.sh

Test a valid file.

  $ cat >valid-export.json <<~~
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
  > }
  > ~~

  $ ${AUGUR} validate export-v2 valid-export.json
  Validating schema of 'valid-export.json'...
  	WARNING:  No colourings were provided
  Validation of 'valid-export.json' succeeded, but there were warnings you may want to resolve.
  Validating that the JSON is internally consistent...
