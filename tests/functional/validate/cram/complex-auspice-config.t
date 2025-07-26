Setup

  $ source "$TESTDIR"/_setup.sh

Test auspice config with invalid colorings:

  $ cat >invalid-colorings.json <<~~
  > {
  >   "title": "Test Dataset",
  >   "colorings": [
  >     {
  >       "key": "none",
  >       "title": "Invalid None Key"
  >     }
  >   ]
  > }
  > ~~

  $ ${AUGUR} validate auspice-config-v2 invalid-colorings.json
  Validating schema of 'invalid-colorings.json'...
    .colorings[0].key failed: 'none' should not be valid under {'const': 'none'}
  FATAL ERROR: Validation of 'invalid-colorings.json' failed.
  [2]

Test auspice config with invalid coloring type:

  $ cat >invalid-coloring-type.json <<~~
  > {
  >   "title": "Test Dataset",
  >   "colorings": [
  >     {
  >       "key": "region",
  >       "title": "Region",
  >       "type": "invalid_type"
  >     }
  >   ]
  > }
  > ~~

  $ ${AUGUR} validate auspice-config-v2 invalid-coloring-type.json
  Validating schema of 'invalid-coloring-type.json'...
    .colorings[0].type failed: 'invalid_type' is not one of ['continuous', 'temporal', 'ordinal', 'categorical', 'boolean']
  FATAL ERROR: Validation of 'invalid-coloring-type.json' failed.
  [2]

Test auspice config with invalid hex color:

  $ cat >invalid-hex-color.json <<~~
  > {
  >   "title": "Test Dataset",
  >   "colorings": [
  >     {
  >       "key": "region",
  >       "title": "Region",
  >       "scale": [
  >         ["Europe", "#invalid"],
  >         ["Asia", "#00FF00"]
  >       ]
  >     }
  >   ]
  > }
  > ~~

  $ ${AUGUR} validate auspice-config-v2 invalid-hex-color.json
  Validating schema of 'invalid-hex-color.json'...
    .colorings[0].scale[0][1] failed: Expected a value with the pattern ^#[0-9A-Fa-f]{6}$ but found "#invalid"
  FATAL ERROR: Validation of 'invalid-hex-color.json' failed.
  [2]

Test auspice config with missing required coloring key:

  $ cat >missing-coloring-key.json <<~~
  > {
  >   "title": "Test Dataset",
  >   "colorings": [
  >     {
  >       "title": "Missing Key"
  >     }
  >   ]
  > }
  > ~~

  $ ${AUGUR} validate auspice-config-v2 missing-coloring-key.json
  Validating schema of 'missing-coloring-key.json'...
    .colorings[0] failed: Missing required field 'key'
  FATAL ERROR: Validation of 'missing-coloring-key.json' failed.
  [2]

Test auspice config with invalid maintainer structure:

  $ cat >invalid-maintainer.json <<~~
  > {
  >   "title": "Test Dataset",
  >   "maintainers": [
  >     {
  >       "name": "Test Author"
  >     }
  >   ]
  > }
  > ~~

  $ ${AUGUR} validate auspice-config-v2 invalid-maintainer.json
  Validating schema of 'invalid-maintainer.json'...

Test auspice config with additional properties:

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
    top level failed: Unexpected property 'invalid_field'
  FATAL ERROR: Validation of 'additional-properties.json' failed.
  [2]
