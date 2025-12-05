Setup

  $ source "$TESTDIR"/_setup.sh

All properties in the schema are optional, so an empty JSON is valid.

  $ cat >minimal-config.json <<~~
  > {}
  > ~~

  $ ${AUGUR} validate auspice-config-v2 minimal-config.json
  Validating schema of 'minimal-config.json'...
