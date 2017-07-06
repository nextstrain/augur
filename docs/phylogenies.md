## building the tree

Building the tree is a 2 step process: A newick file is created (e.g. with RAxML) and then this is converted into a TimeTree.

* `tree_verbosity_level` {int} (default `0`). Control the amount of information displayed. Higher numbers -> more output.

### step 1: initial phylogeny
This is normally created via RAxML however FastTree may also be used. A number of options are available however for the most part the defaults are fine.
All the keyword arguments in the `build_newick` method of the `Tree` class may be set via the `newick_tree_options` config dictionary. The most common ones are:

* `nthreads` {int} (default `2`) how many threads will RAxML use
* `num_distinct_starting_trees` {int} (default `1`) more is better but sloooooower

By default **RAxML** is run using the GTRCAT model with a single starting tree.

The newick file is saved to `processed/<prefix>.newick`. If this file exists, and the taxa match, then this step is skipped to save time.

### step 2: TimeTree
TimeTree involves three steps:
* set up
* clock filtering (configurable via `config -> clock_filter`)
* reconstruct ancestral sequences, date nodes & optimise branch lengths (configurable via `config -> timetree_options`)

**clock_filter:**

`clock_filter` is a {dict} with keys:
* `n_iqd` {int} (default `3`)
* `plot` {bool} (default `True`)
* `remove_deep_splits` {bool} (default `False`)

Note that setting `clock_filter` to `False` skips this step.

**timetree_options:**
Available parameters:
* `Tc` (default `0.02`)
* `confidence` (see below)
* `use marginal` (see below)
* `reroot` (default `"best"`)
* `resolve_polytomies`
* others related to skyline (todo)

Note that setting `config["temporal confidence"]` turns on `confidence` and `use_marginal`
