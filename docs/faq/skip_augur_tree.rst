=================================
How do I use my own tree builder?
=================================

The `augur tree` command is a light wrapper around tree building programs such as IQ-TREE, RAxML and FastTree.
It's possible that the functionality you want isn't available in those programs, or that it is available but that `augur tree` doesn't expose the functionality you need.
The solution is just to use your own tree building program directly.
This guide explains the general approach you can use.


Exclude sites
-------------

Augur tree allows you to mask certain problematic sites before tree-building via  `augur tree --exclude-sites <exclude_file> --alignment <alignment>` , where each line in `<exclude_file>` is a 1-based nucleotide position to mask from the alignment. We can recreate this via the following command:

.. code-block:: bash

  augur mask --sequences <alignment> --mask <exclude_file> --output <masked_alignment>


If you use this, remember to use the `<masked_alignment>` file for tree building, not the input `<alignment>` file!

Build your tree
---------------

This is up to you!
The end result should be a newick tree file with branch lengths representing nucleotide divergence (or similar).
In a typical nextstrain workflow, this file will be handed to `augur refine` for temporal inference (see below if you also want to skip this step).
Note that augur refine may change the tree structure in two ways:

* Rerooting (avoid via `augur refine --keep-root`)
* Resolving polytomies using temporal information  (avoid via `augur refine --keep-polytomies`)


There are some gotchas involved in tree building, such as:

* IQ-TREE will change taxon names if they include any of the following characters: `/|()*`. This will result in the tree and alignment falling out of sync. Other tree builders may have similar behavior.
* Bootstrap values are allowed, and should be parsed correctly by `augur refine`.



Avoiding augur refine
---------------------

While augur refine is typically used to reroot the tree and infer the times of internal nodes, it also performs two important tasks which subsequent augur commands expect:

* Labels internal node names. These are typically in the format `NODE_0000001`, but this is not required.
* Stores branch lengths already present in the newick tree as a node-data JSON file, which `augur export` uses.

If you can address both of these then you can easily skip the `augur refine` step.

