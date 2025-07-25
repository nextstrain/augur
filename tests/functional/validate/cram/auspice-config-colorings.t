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
    .colorings[0].key "none" failed not validation for {"const": "none"}
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
    .colorings[0].type "invalid_type" failed enum validation for ["continuous", "temporal", "ordinal", "categorical", "boolean"]
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
    .colorings[0].scale[0][1] "#invalid" failed pattern validation for "^#[0-9A-Fa-f]{6}$"
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
    .colorings[0] {"title": "Missing Key"} failed required validation for ["key"]
  FATAL ERROR: Validation of 'missing-coloring-key.json' failed.
  [2]
