Setup

  $ source "$TESTDIR"/_setup.sh

Test a valid maintainers entry.

  $ cat >valid-auspice-config.json <<~~
  > {
  >   "title": "Test Dataset",
  >   "maintainers": [
  >     {"name": "Test Author", "url": "https://example.com"}
  >   ]
  > }
  > ~~

  $ ${AUGUR} validate auspice-config-v2 valid-auspice-config.json
  Validating schema of 'valid-auspice-config.json'...
