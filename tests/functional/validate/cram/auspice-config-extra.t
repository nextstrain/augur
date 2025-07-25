Setup

  $ source "$TESTDIR"/_setup.sh

Additional properties are not allowed.

  $ cat >additional-properties.json <<~~
  > {
  >   "title": "Test Dataset",
  >   "invalid_field": "should not be here",
  >   "maintainers": [
  >     {"name": "Test Author", "url": "https://example.com"}
  >   ]
  > }
  > ~~

  $ ${AUGUR} validate auspice-config-v2 additional-properties.json
  Validating schema of 'additional-properties.json'...
     {"title": "Test Dataset", "invalid_field": "shou…} failed additionalProperties validation for false
  FATAL ERROR: Validation of 'additional-properties.json' failed.
  [2]
