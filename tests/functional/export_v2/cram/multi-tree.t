
Setup

  $ source "$TESTDIR"/_setup.sh

Create a small second tree (which has different names/labels than 'tree.nwk')
  $ cat > tree2.nwk <<~~
  > (tipG:1,(tipH:1,tipI:1)internalHI:2)SECOND_ROOT:0;
  > ~~

Minimal export -- no node data, no metadata etc etc
  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" tree2.nwk \
  >   --output minimal.json &> /dev/null

More realistic export - with node_data for all nodes and metadata for some of them

  $ cat > metadata.tsv <<~~
  > strain	something
  > tipA	foo
  > tipB	foo
  > tipC	foo
  > tipG	bar
  > tipH	bar
  > ~~


  $ cat > node-data.json <<~~
  > {
  >   "nodes": {
  >     "ROOT": {"mutation_length": 0},
  >     "tipA": {"mutation_length": 1},
  >     "internalBC": {"mutation_length": 2},
  >     "tipB": {"mutation_length": 1},
  >     "tipC": {"mutation_length": 1},
  >     "internalDEF": {"mutation_length": 5},
  >     "tipD": {"mutation_length": 3},
  >     "tipE": {"mutation_length": 4},
  >     "tipF": {"mutation_length": 1},
  >     "SECOND_ROOT": {"mutation_length": 0},
  >     "tipG": {"mutation_length": 1},
  >     "internalHI": {"mutation_length": 2},
  >     "tipH": {"mutation_length": 1},
  >     "tipI": {"mutation_length": 1}
  >   }
  > }
  > ~~

  $ ${AUGUR} export v2 \
  >   --tree "$TESTDIR/../data/tree.nwk" tree2.nwk \
  >   --metadata metadata.tsv --color-by-metadata something \
  >   --node-data node-data.json \
  >   --output output.json &> /dev/null
