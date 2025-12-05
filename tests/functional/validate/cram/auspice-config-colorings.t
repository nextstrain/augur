Setup

  $ source "$TESTDIR"/_setup.sh

Test a valid colorings entry.

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

Test with invalid colorings.

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

Test with invalid coloring type.

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

Test with invalid hex color.

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

Test with missing required coloring key.

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
    .colorings[0] failed: Missing required property 'key'
  FATAL ERROR: Validation of 'missing-coloring-key.json' failed.
  [2]
