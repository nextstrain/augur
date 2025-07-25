Setup

  $ source "$TESTDIR"/_setup.sh

Test valid auspice config:

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

Test minimal valid auspice config:

  $ cat >minimal-config.json <<~~
  > {}
  > ~~

  $ ${AUGUR} validate auspice-config-v2 minimal-config.json
  Validating schema of 'minimal-config.json'...

Test valid auspice config with colorings:

  $ cat >config-with-colorings.json <<~~
  > {
  >   "title": "Test Dataset with Colorings",
  >   "colorings": [
  >     {
  >       "key": "region",
  >       "title": "Region",
  >       "type": "categorical"
  >     }
  >   ]
  > }
  > ~~

  $ ${AUGUR} validate auspice-config-v2 config-with-colorings.json
  Validating schema of 'config-with-colorings.json'...
