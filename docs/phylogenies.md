## building the tree

Building the tree is a 2 step process: A newick file is created (e.g. with RAxML) and then this is converted into a TimeTree.

### step 1: initial phylogeny
This is normally created via RAxML however FastTree may also be used. A number of options are available however for the most part the defaults are fine.
All the keyword arguments in the `build_newick` method of the `Tree` class may be set via the `newick_tree_options` config dictionary.

By default **RAxML** is run using the GTRCAT model with a single starting tree.

The newick file is saved to `processed/<prefix>.newick`. If this file exists, and the taxa match, then this step is skipped to save time.

### step 2: TimeTree
Options here are configurable via the `timetree_options` dictionary in config.

Available parameters:
* `Tc` (default `0.02`)
* `confidence` (see below)
* `use marginal` (see below)
* `reroot` (default `"best"`)
* `resolve_polytomies`
* others related to skyline (todo)

Note that setting `config["temporal confidence"]` turns on `confidence` and `use_marginal`
