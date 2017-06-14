## Augur pipeline steps

## Initialization
Sets up the main object with a number of (modifiable) parameters such as paths to the data. If possible, the reference sequence is loaded (using BioPython) and the proteins are identified as features on this reference. Also creates the data/build/store folders in the directory where you run the script.

#### call:
`zika = process()`

#### inputs (all optional):
* `input_data_path` [default: 'data/test'] the file 'data/test.fasta' must exist (this is the main input!)
* `store_data_path` [default: 'store/test']
* `build_data_path` [default: 'build/test']
* `verbose` [default: 2]
* any number of further keyword arguments such as:
  * `reference`: path to reference file
  * `lat_long_fname`: path to lat/long tsv file
  * `proteins`: list of protein names (found in reference)

#### returns:
an object of class `Process`

#### side-effects:
* The following parameters are stored as top level parameters (i.e. `zika.x`): `input_data_path`, `store_data_path`, `build_data_path`, `verbose`
* All other parameters are stored in `zika.kwargs`
* If the input/store/build folders don't yet exist, they are created
* `zika.lat_long_fname` is set if provided as an option, and is necessary for geo-inference
* By default `zika.proteins={}` and `zika.reference_seq=None`
* If a reference is provided then it's loaded as a 'Bio.SeqRecord.SeqRecord' object at `zika.reference_seq`. The `id` of the object is set to its `name`. `zika.genome_annotation` is the the feature list of the reference.
* If a reference and list of proteins are provided, `zika.proteins` is a map of gene names to FeatureLocation objects.
* `zika.sequence_fname` is set to the `input_data_path` + `.fasta`
* `zika.file_dumps` is a map of 3 keys `seqs`, `tree`, `nodes` and paths (in the store directory) to (not yet created) filenames where progress will be saved.








## load the input data
The data (fasta) is loaded into an object from the `sequence_set` class. The sequences are ungapped and the dates are parsed. By default sequences without dates are removed.

#### call:
`zika.load_sequences()` (_process_ method)

#### inputs:
* `fields`: a map of fasta header fields and the keyword for them to be mapped to. [Default: `{0:'strain', 2:'isolate_id', 3:'date', 4:'region',                     5:'country', 7:"city", 12:"subtype",13:'lineage'}`]
* `prune`: should sequences without dates be removed? [Default: True]
* `sep`: fasta seperator character [default: '|']

#### returns:
None

#### side-effects:
* `zika.seqs` is an object from _sequence_set_
  * `zika.seqs.all_seqs` is a map from _seq.name_ -> _SeqRecord_ for all sequences
  * `zika.seqs.run_dir` is a temporary directory created
  * `zika.seqs.reference_seq` is a link to the _SeqRecord_ of the reference (if it's in the fasta the link will be to `zika.seqs.all_seqs[X]`, else to `zika.reference_seq`), else it's `None`
  * all sequences (in `zika.seqs.all_seqs`) are ungapped and capatilised
  * for each sequence, the desired fields are set (e.g. `zika.seqs.all_seqs['ABC'].country="Brazil"`). If _strain_ is one of these fields, then `.id` and `.name` are set to this. If a field is not found then the value is an empty string.
* assuming `date` is a field and the format is `"%Y-%m-%d"`, these attributes on each sequence are created / changed:
  * `.date`: the string from the FASTA file
  * `.num_date`: float e.g. 2016.53, or `[min, max]` if ambiguity
  * `.raw_date`: a _datetime_ object (or, if this fails, the original string)
* all sequences (`zika.seqs.all_seqs`) without a corectly parsed date are removed
* some STDOUT


#### notes:
A number of *sequences_set* methods accomplish this.










## Filter sequences
Sequences are filtered here based upon the user supplying a function which returns a boolean determining whether to keep the sequence or not! This is usually done multiple times to filter on multiple criteria!

#### call:
zika.seqs.filter() (*sequences_set* method)

#### inputs:
* `func`: function returning a truthy value. Its input is a _SeqRecord_ object.
* `leave_ref`: should the reference be forced to stay? [default: False]

#### returns:
None

#### side-effects:
* sequences are pruned from the `zika.seqs.all_seqs` dict.
* some STDOUT

#### notes:
`zika.py` has the following filtering steps: Date range is specified, sequences must be over 2kb, and a list of strains to drop is specified.







## Subsample Sequences
This uses user-defined functions to restrict the data to a manageable size. Based upon these functions, sequences are assigned to categories for subsampling, assigned priorities within each category, sorted upon this priority and then the threshold determines how many to take. _Note that the sequences are now stored in a different location!_

#### call:
zika.seqs.subsample() (*sequences_set* method)

#### inputs:
* `category`: callable (arg: _SeqRecord_ object, can return pretty much anything to use as an ID)
* `priority`: callable (arg: _SeqRecord_)[default: random priority]
* `threshold`: INT or callable (args: the category ID and the _SeqRecord_) [default: 5]
* `repeated`: set this if you're running this function subsequent times [default: False]
* `forced_strains`: list of of strain names that should always be included [default: empty list]

#### returns:
None

#### side-effects:
* `zika.seqs.seqs` is the subsampled version of `zika.seqs.all_seqs`, that is, a map of sequence name to _SeqRecord_ (links)
* `zika.seqs.sequence_categories` is a dict of category ID (see inputs) to a list of _SeqRecord_ objects (these are links)
* each _SeqRecord_ has the attribute `under_sampling` (FLOAT) set, which is actually the degree of over sampling in each category as well as `_priority` (FLOAT)

#### notes:
To randomly subsample 5 viruses per month, you would do something like `zika.seqs.subsample(category = lambda x:(x.attributes['date'].year, x.attributes['date'].month), threshold=5)`












## Align and Translate Sequences
Alignments use _mafft_ and include the reference sequence (if it wasn't already there). For codon alignments, sequences with premature stops are removed. Positions with gaps in the reference alignment are removed from all the alignments.

If an outgroup is set, remove sequences form the set that are that evolve much faster or slower compared the majority. Regions with predominantly gaps can be removed since this can skew the evolutionary rates. (The reference seq will not be removed.)

Finally, each protein sequence is translated (if there were no proteins set, the whole alignment is translated).

#### call:
`zika.align()` (_Process_ method)

#### inputs:
* `outgroup` [Default: None]
* `codon_align` [Default: False]
* `debug` keeps the temporary folder with the alignmetns [Default: False]

#### returns:
None

#### side-effects:
* creates the `run_dir` and (unless debug is True) deletes it
* `zika.seqs.aln` the nucleotide alignment (with `-` as gaps) (a _MultipleSeqAlignment_ object). Attributes from `zika.seqs.seqs` are linked here.
* `zika.seqs.sequence_lookup` a map linking seq.id to the alignment in `zika.seqs.aln`
* a matplotlib plot is created if outgroup is set.
* `zika.seqs.translations` is a dict with keys of the proteins (cds if none set) and values objects of _Bio.Align.MultipleSeqAlignment_
* `zika.seqs.proteins` is created from `zika.proteins`



#### notes:
This will throw if the reference is not set.











## Build Tree
FastTree is used to create an initial tree, then RAxML refines the topology (under a GTRCAT model) and then refines the branch lengths (GTRGAMMA model). This tree is then handed to the _treetime_ module, which re-roots the tree, finds clades and (if set) adds node properties from a provided file.


#### call:
`zika.build_tree()`

#### inputs:
* infile: If provided, skips the tree building steps. [Default: None]
* nodefile a cPickle file containing a dict of node name -> dict of attributes which will be added to the tree [Default: None]
* root [Default: 'best']
* debug: controls whether temporary files should be removed [Default: False]

#### returns:
None

#### side-effects:
* `zika.tree` is a _tree_ object, with the following important properties:
  * `zika.tree.aln` is a link to `zika.seqs.aln`
  * `zika.tree.proteins` is a link to `zika.proteins`
* Temporary directory is created (in the directory the script was called from)
* `zika.tree.tt` is a _treetime_ object, and the tree is at `zika.tree.tt.tree`
* `zika.tree.tree` is a link to `zika.tree.tt.tree`
* `zika.tree.is_timetree = False`
* The tree has the following node properties:
  * the alignment sequences added to the leaves.
  * Dates are added to each leaf (`node.numdate_given` and `node.date`.
  * If a date is not available then `node.bad_branch` is `True`. Internal nodes where all leaves don't have dates have this property set as well.
  * Any properties in the provided _nodefile_ are added
* `zika.tree.dump_attr` is an empty list


#### notes:
FastTree and RAxML binaries must be installed and callable via `fasttree` and `raxml`. The tree itself is a [Phylo](http://biopython.org/DIST/docs/tutorial/Tutorial.html) object, with the important method `find_clades` which searches the tree and returns an iterable through all matching objects.


Requires [treetime](https://github.com/neherlab/treetime). If you're reading that code you might want to take a look at `property` decoratation in python ([see here](http://stackoverflow.com/questions/2627002/whats-the-pythonic-way-to-use-getters-and-setters)). `zika.tree.tt` is an object of class _TreeTime_, which inherits from _ClockTree_ which inherits from _TreeAnc_.








## Save Progress
The attributes specified in `zika.file_dumps` (the keys) are saved in the filenames (the values of the dict) using cPickles. The `nodes` attribute is created (actually called `node_props`) as a dictionary by crawling the tree.

#### call:
`zika.dump()`

#### inputs:
None

#### returns:
None

#### side-effects:
File creation










## Molecular Clock Filtering
This method labels outlier branches that don't seem to follow a molecular clock and excludes them from subsequent molecular clock estimate and timetree propagation.


#### call:
`zika.clock_filter`

#### inputs:
* `n_iqd` [default: 3]
* `plot` [default: True]

#### returns:
None

#### side-effects:
* If `plot` is truthy, a root-to-tip plot is generated.
* Leaves are pruned (according to the `node.bad_branch` boolean which is set during this method, so only nodes with Truthy values remain)
* The tree is re-rooted, which results in topology and attribute changes to `zika.tree.tt` (such as `.tree.root.branch_length`, `tree.root.mutations`, internal clade names on the tree etc)











## Annotate Tree
Now that we have a tree, this method re-roots the tree, optimizes the branch lengths and reconstructs ancestral sequences at each node. Optionally, the nodes are dated, including optional confidence estimation via marginal reconstruction. The seqeunces at each node are translated (using the protein data supplied) and AA mutations are found. To prepare for export, the tree is ladderized.

#### call:
`zika.annotate_tree()`

#### inputs:
* `Tc` The coalescent prior. Can be float, 'opt' or 'skyline' [default: 0.01]
* `timetree`: bool. Run node dating? [Default: False]
* `kwargs` including:
  * `reroot` [Default: 'best']
  * `resolve_polytomies` [Default: True]
  * `max_iter` Number of timetree iterations [Default: 2]
  * `confidence` Run dating confidence estimation [Default: False]


#### returns:
None

#### side-effects:
These methods cause a huge number of side-effects within the timetree object but very few elsewhere.
* `zika.tree.is_timetree` is set to true if node dating was run
* `zika.tree.dump_attr` list has the following strings added:
  * `sequence` (always)
  * `numdate`, `date` (if node dating was run)
  * `translations` (always)
  * `muts`, `aa_muts`, `aa_mutations`, `mutation_length`, `mutations`
  * `yvalue`, `xvalue`, `clade` (always)
  * `tvalue` (if node dating was run)
* each `node` on the tree (`zika.tree.tt.tree` or `zika.tree.tree`) now has:
  * `node.sequence` (via `optimize_sequences_and_branch_length`)
  * `node.mutations` a list of tuples (a,pos,d), where a is the parental nt and d the nt at this node. Pos is 0-based. (via `optimize_sequences_and_branch_length`)
  * `node.attr` has `num_date` and (optionally) `num_date_confidence` set (via the `timetree` method)
  * `node.translations` a dict of protein name to (inferred) protein sequence (a _SeqFeature_ object) (via `add_translations`)
  * `node.aa_mutations` is a dictionary of protein names to a list of tuples similar to `node.mutations` but referring to AA not nt (via `add_translations`)
  * `node.muts` is a string representation of `node.mutations` but is 1-based not 0-based (via `refine` method of _tree_)
  * `node.aa_muts` is, similarly, the string representation of the tuples in `node.aa_mutations` (via `refine` method of _tree_)
  * `node.attr['div']` is the branch length _from the root_ (via `refine` method of _tree_)
  * `node.xvalue` the same as `node.attr['div']` (via `layout` method of _tree_)
  * `node.tvalue` time value since the root (via `layout` method of _tree_)
  * `node.yvalue` For leaves, this is an integer referring to tip ordering, using the "preorder" method of `find_clades` for the leaves, for internal nodes its the mean value of all of the children tips (via `layout` method of _tree_)



#### notes:
* `node.up` must refer to the parent (under the "preorder" method), as leaf nodes still have this
* Nodes on the tree have the attribute `attr` set, which is a dictionary. This is confusing as, e.g. if `node.attr.region = 'Brazil'`, `hasattr(node, 'region')` is False.
* Inspect a (leaf) node's attributes via `zika.tree.tree.get_terminals()[0].attr`
* The `dump_attr` list controls what's saved if `zika.dump()` is run.





## Geographical Inference
Infers a "mugration" model by pretending each region corresponds to a sequence state and repurposing the GTR inference and ancestral reconstruction. Normally this method is called in a loop over a list.

#### call:
`zika.tree.geo_inference()`

#### inputs:
A string which should be found in `node.attr` for the leaves of the tree. E.g. "region" or "country".

#### returns:
None

#### side-effects:
A number of attributes of `zika.tree.tt` are modified. Namely:
* The input string is added to `zika.tree.dump_attr`
* For each `node` in the tree:
  * The input string is a key in `node.attr` with the inferred value
  * `node.X_transitions` (where `X` is the input string) is the geo-transitions, in the same format as `node.mutations`

#### notes:
* This involves temporarily re-defining (i.e. changing) the properties (such as sequence, mutations) of the nodes on the tree. They are set back before this method returns.




## Export Data
This methods has a rather complicated set of inputs which control the data exported for auspice.

#### call:
`zika.export`

#### inputs:
* `extra_attr` list of extra attributes to save with the tree [Default: None]
* `controls` A dict of "output name" and lists of "fields" are used to control the meta JSON output. The "output name" is the JSON field (within "controls"), and the number of terminal nodes in the tree with `field in node.attr` are reported. E.g. `{'geographic location':['region', 'country']}`. [Default: empty dict].
* `geo_attributes` Should be an array of strings for which `geo_inference` was run. [Default: empty list]
* `color_options` Saved directly to meta JSON output [Default: `{"num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"}}`]
* `panels` Saved directly in meta JSON output [Default = `['tree', 'entropy']`]
* `indent` Used for tree output formatting [Default: None]

#### returns:
None

#### side-effects:
A number of files are created (as you would expect!) and saved in the `build` directory (see `zika.build_data_path`)
* Alignment diversity (`entropy.json`) is a JSON file that contains alignment diversity column by column
* Tree (`tree.json`) with attributes `["muts", "aa_muts","attr", "clade"]`
* Frequency estimations (`frequencies.json`) if `zika.estimate_mutation_frequencies` was run, that is if `zika.pivots` exists. The format is region_protein:159F for mutations and region_clade:123 for clades.
* Metadata (`meta.json`) including some input data (see above), the git commit, today's date, the controls (see above), and geo attributes which is a mapping from the countries (etc) to lat long points.

#### notes:
* This requires the file set at `zika.lat_long_fname` to exist, which was defined when the `zika` object was created. This would be best checked for when `zika` is created.
* The only (currently) exportable strings from `geo_attributes` are "region", "country" and "division".
