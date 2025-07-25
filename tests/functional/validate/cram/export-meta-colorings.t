Setup

  $ source "$TESTDIR"/_setup.sh

Test export with no colorings.

  $ cat >no-colorings.json <<~~
  > {
  >   "version": "v2",
  >   "meta": {
  >     "updated": "2023-12-25",
  >     "panels": ["tree"]
  >   },
  >   "tree": {
  >     "name": "root",
  >     "node_attrs": {"div": 0, "num_date": {"value": 2023.0} }
  >   }
  > }
  > ~~

  $ ${AUGUR} validate export-v2 no-colorings.json
  Validating schema of 'no-colorings.json'...
  	WARNING:  No colourings were provided
  Validation of 'no-colorings.json' succeeded, but there were warnings you may want to resolve.
  Validating that the JSON is internally consistent...

Test with missing colorings reference in meta.

  $ cat >missing-coloring-ref.json <<~~
  > {
  >   "version": "v2",
  >   "meta": {
  >     "updated": "2023-12-25",
  >     "panels": ["tree"],
  >     "colorings": [
  >       {
  >         "key": "region",
  >         "title": "Region"
  >       }
  >     ]
  >   },
  >   "tree": {
  >     "name": "root",
  >     "node_attrs": {"div": 0, "num_date": {"value": 2023.0} },
  >     "children": [
  >       {
  >         "name": "child1",
  >         "node_attrs": {"div": 0.1, "num_date": {"value": 2023.1} }
  >       }
  >     ]
  >   }
  > }
  > ~~

  $ ${AUGUR} validate export-v2 missing-coloring-ref.json
  Validating schema of 'missing-coloring-ref.json'...
    .meta.colorings[0] {"key": "region", "title": "Region"} failed required validation for ["key", "type"]
  FATAL ERROR: Validation of 'missing-coloring-ref.json' failed.
  [2]
