Setup

  $ source "$TESTDIR"/_setup.sh

Test duplicate node names.

  $ cat >duplicate-nodes.json <<~~
  > {
  >   "version": "v2",
  >   "meta": {
  >     "updated": "2023-12-25",
  >     "panels": ["tree"]
  >   },
  >   "tree": {
  >     "name": "root",
  >     "node_attrs": {"div": 0, "num_date": {"value": 2023.0} },
  >     "children": [
  >       {
  >         "name": "SEQ1",
  >         "node_attrs": {"div": 0.1, "num_date": {"value": 2023.1} }
  >       },
  >       {
  >         "name": "SEQ1",
  >         "node_attrs": {"div": 0.2, "num_date": {"value": 2023.2} }
  >       }
  >     ]
  >   }
  > }
  > ~~

  $ ${AUGUR} validate export-v2 duplicate-nodes.json
  Validating schema of 'duplicate-nodes.json'...
  FATAL ERROR: Node SEQ1 appears multiple times in the tree.
  Validating that the JSON is internally consistent...
  [2]
