# CHANGELOG

## __NEXT__

### Bug Fixes

* mask: Fix parsing of BED files as zero-indexed, half-open intervals [#512][]

### Features

* mask: Allow masking of specific sites passed by the user with `--mask-sites` and masking of a fixed number of sites from the beginning or end of each sequence with `--mask-from-beginning` and `--mask-from-end` [#512][]

[#512]: https://github.com/nextstrain/augur/pull/512

## 7.0.2 (7 April 2020)

### Bug Fixes

* filter: Fix regression introduced in 7.0.0 which caused an error to be raised
  if a priorities file didn't include every sequence.  Sequences which are not
  explicitly listed will once again default to a priority of 0. [#530][]

[#530]: https://github.com/nextstrain/augur/pull/530


## 7.0.1 (7 April 2020)

### Bug Fixes

* Fix typo with Python classifiers in setup.py

## 7.0.0 (7 April 2020)

### Major Changes

* Drop support for Python 3.4 and 3.5 [#482][]
* Drop support for `--output` flag in augur ancestral, clades, sequence-traits, traits, and translate in favor of `--output-node-data` flag [#529][]

### Features

* improve testing by
	* adding a simple shell script to run tests and improving pytest configuration and output [#463][]
	* adding code coverage reports ([#486][], [#491][]) and integration with codecov.io [#508][]
	* adding unit tests for align ([#477][]), filter ([#478][], [#479][], [#487][]), utils ([#501][])
* align: reverse complement sequences when necessary using mafft’s autodirection flag [#467][]
* align: speed up replacement of gaps with “ambiguous” bases [#474][]
* mask: add support for FASTA input files [#493][]
* traits: bump TreeTime version to 0.7.4 and increase maximum number of unique traits allowed from 180 to 300 [#495][]

### Bug Fixes

* align: enable filling gaps in input sequences even if no reference is provided instead of throwing an exception [#466][]
* align: detect duplicate sequences by comparing sequence objects instead of (often truncated) string representations of those objects [#468][]
* import_beast: use raw strings for regular expressions to avoid syntax errors in future versions of Python [#469][]
* scripts: update exception syntax to new style [#484][]
* filter: fail loudly when a given priority file is invalid and exit instead of just printing an error [#487][]

### Documentation

* README: document need for brew tap prior to brew install [#471][]
* DEV_DOCS: add a proper section on how to write and run tests [#481][]
* README: improve discoverability of core documentation with RTD badge [#490][]

[#463]: https://github.com/nextstrain/augur/pull/463/
[#466]: https://github.com/nextstrain/augur/pull/466
[#467]: https://github.com/nextstrain/augur/pull/467
[#468]: https://github.com/nextstrain/augur/pull/468
[#469]: https://github.com/nextstrain/augur/pull/469
[#471]: https://github.com/nextstrain/augur/pull/471
[#474]: https://github.com/nextstrain/augur/pull/474
[#477]: https://github.com/nextstrain/augur/pull/477
[#478]: https://github.com/nextstrain/augur/pull/478
[#479]: https://github.com/nextstrain/augur/pull/479
[#481]: https://github.com/nextstrain/augur/pull/481
[#482]: https://github.com/nextstrain/augur/pull/482
[#484]: https://github.com/nextstrain/augur/pull/484
[#486]: https://github.com/nextstrain/augur/pull/486
[#487]: https://github.com/nextstrain/augur/pull/487
[#490]: https://github.com/nextstrain/augur/pull/490
[#491]: https://github.com/nextstrain/augur/pull/491
[#493]: https://github.com/nextstrain/augur/pull/493
[#495]: https://github.com/nextstrain/augur/pull/495
[#501]: https://github.com/nextstrain/augur/pull/501
[#508]: https://github.com/nextstrain/augur/pull/508
[#529]: https://github.com/nextstrain/augur/pull/529

## 6.4.3 (25 March 2020)

### Bug Fixes

* align: Remove reference sequence from alignments even when no gaps exist in
  input sequences relative to the reference. Thank you @danielsoneg!
  [#456][]

### Documentation

* Reorganize README, improve findability of documentation, and add separate dev
  docs.
  [#461][]

[#456]: https://github.com/nextstrain/augur/pull/456
[#461]: https://github.com/nextstrain/augur/pull/461

## 6.4.2 (17 March 2020)

### Bug Fixes

* Require Snakemake less than 5.11 to avoid a breaking change.  The `--cores`
  argument is now required by 5.11, which will affect many existing augur-based
  workflows.  Reported upstream as
  [snakemake/snakemake#283](https://github.com/snakemake/snakemake/issues/283).

* align: Run mafft with the `--nomemsave` option.  This makes alignments of
  sequences over 10k in length run much, much faster in the general case and
  shouldn't cause issues for most modern hardware.  We may end up needing to
  add an off-switch for this mode if it causes issues for other users of augur,
  but the hope is that it will make things just magically run faster for most
  folks!  There is likely more tuning that could be done with mafft, but this
  is a huge improvement in our testing.
  [#458][]

* align: Ignore blank lines in `--include` files.  Thanks @CameronDevine!
  [#451][]

* align: Properly quote filenames when invoking mafft.  Thanks @CameronDevine!
  [#452][]


[#451]: https://github.com/nextstrain/augur/pull/451
[#452]: https://github.com/nextstrain/augur/pull/452
[#458]: https://github.com/nextstrain/augur/pull/458


## 6.4.1 (4 March 2020)

### Bug Fixes

* export: AA labels are now exported for branches where a clade is also labeled
  [See PR 447](https://github.com/nextstrain/augur/pull/447)

* export / validation: a dataset title is no longer required

* release script now works on MacOS & code-signing is optional
  [See PR 448](https://github.com/nextstrain/augur/pull/448)

* traits: Missing data is correctly handled

## 6.4.0 (26 February 2020)

### Features

* align: New sequences can now be added to an existing alignment.  [#422][]

* align: Multiple sequence files can be provided as input. [#422][]

* align: Extra debugging files such as `*.pre_aligner.fasta` and
  `*.post_aligner.fasta` are no longer produced by default.  To request them,
  pass the `--debug` flag. [#422][]

* align: De-duplicate input sequences, with a warning. [#422][]

* export v2: Add support for the `branch_label` property in `display_defaults`,
  which was recently added to Auspice. [#445][]

### Bug fixes

* align: Exits with an error earlier if arguments are invalid instead of only
  printing a warning. [#422][]

* align: Performs more error checking and clarifies the help and error
  messages. [#422][]

* export v2: Traits which are filters but not colorings are now exported as
  well, instead of being left out. [#442][]

* export v2: Exits non-zero when validation fails, instead of masking errors.
  [#441][]

* validate: In order to improve clarity, messages now include the filenames
  involved and distinguish between schema validation and internal consistency
  checks. [#441][]

[#422]: https://github.com/nextstrain/augur/pull/422
[#441]: https://github.com/nextstrain/augur/pull/441
[#442]: https://github.com/nextstrain/augur/pull/442
[#445]: https://github.com/nextstrain/augur/pull/445


## 6.3.0 (13 February 2020)

### Features

* Augur `refine`, `ancestral` and `traits` now use the
  [upgraded TreeTime v0.7](https://github.com/neherlab/treetime/releases/tag/v0.7.0)
  This should have a number of under-the-hood improvements.
  [See PR 431](https://github.com/nextstrain/augur/pull/431)
* ancestral: New options to either `--keep-ambiguous` or `--infer-ambiguous`. If using
  `--infer-ambiguous` the previous behavior will be maintained in which tips with `N` will have
  their nucleotide state inferred. If using `--keep-ambiguous`, these tips will be left as `N`.
  With this upgrade, we are still defaulting to `--infer-ambiguous`, however, we plan to swap
  default to `--keep-ambiguous` in the future. If this distintion matters to you, we would suggest
  that you explicitly record `--keep-ambiguous` / `--infer-ambiguous` in your build process.
  [Also part of PR 431](https://github.com/nextstrain/augur/pull/431)
* traits: Allow input of `--weights` which references a `.tsv` file in the following format:
  ```
  division	Hubei	10.0
  division	Jiangxi	1.0
  division	Chongqing	1.0
  ```
  where these weights represent equilibrium frequencies in the CTMC transition model. We imagine the
  primary use of user-specified weights to correct for strong sampling biases in available data.
  [See PR 443](https://github.com/nextstrain/augur/pull/443)

### Bug fixes

* Improvements to make shell scripts run more easily on Windows.
  [See PR 437](https://github.com/nextstrain/augur/pull/437)

## 6.2.0 (25 January 2020)

### Features

* refine: Include `--divergence-units` option to distinguish between `mutations`
  and `mutations-per-site`. Keep `mutations-per-site` as default behavior.
  [See PR 435](https://github.com/nextstrain/augur/pull/435)

### Bug fixes

* utils: Support v2 auspice JSONs in json_to_tree utility function.
  [See PR 432](https://github.com/nextstrain/augur/pull/432)

## 6.1.1 (17 December 2019)

### Bug fixes

* frequencies: Fix bug in string matching for weighted frequencies introduced in
  v6.1.0. [See PR 426](https://github.com/nextstrain/augur/pull/426).

## 6.1.0 (13 December 2019)

### Features

* export: Include `--description` option to pass in a Markdown file with dataset
  description. This is displays in Auspice in the footer. For rationale,
  [see Auspice issue 707](https://github.com/nextstrain/auspice/issues/707) and
  for Augur changes [see PR 423](https://github.com/nextstrain/augur/pull/423).

### Bug fixes

* frequencies: Fix weighted frequencies when weight keys are unrepresented.
  [See PR 420](https://github.com/nextstrain/augur/pull/420).


## 6.0.0 (10 December 2019)

### Overview

Version 6 is a major release of augur affecting many augur commands. The format
of the exported JSON (v2) has changed and now merges the previously separate
files containing tree and meta information. To maintain backward compatibility,
the export command was split into `export v1` (old) and `export v2` (new).
Detailed release notes are provided in the augur documentation [on
read-the-docs](https://nextstrain-augur.readthedocs.io/en/stable/releases/v6.html).
For a migration guide, consult
[migrating-v5-v6](https://nextstrain-augur.readthedocs.io/en/stable/releases/migrating-v5-v6.html).

### Major features / changes

* export: Swap from a separate `_tree.json` and `_meta.json` to a single
  "unified" `dataset.json` output file
* export: Include additional command line options to alleviate need for Auspice
  config
* export: Include option for reference sequence output
* export: Move to GFF-style annotations
* export: Validate exported JSONs against schema
* ancestral: Allow output of FASTA and JSON files
* import: Include `import beast` command to import labeled BEAST MCC tree
* parse: Include `--prettify-fields` option to cleanup metadata fields
* Documentation improvements

### Minor features / changes

* colors.tsv: Allow whitespace, but insist on tab delimiting
* lat_longs.tsv: Allow whitespace, but insist on tab delimiting
* Remove code for old "non-modular" augur, old "non-modular" builds and Python
  tests
* Improve test builds
* filter: More interpretable output of how many sequences have been filtered
* filter: Additional flag `--subsample-seed` to seed the random number generator
  and thereby make subsampling reproducible
* sequence-traits: Numerical output as originally intended, but required an
  Auspice bugfix
* traits: Explanation of what is considered missing data & how it is interpreted
* traits: GTR models are exported in the output JSON for better accountability &
  reproducibility


## 5.4.1 (12 November 2019)

### Bug fixes

* export v1: Include `--minify-json` option that was mistakenly not included in PR 398.
  [See PR 409](https://github.com/nextstrain/augur/pull/409)

## 5.4.0 (7 November 2019)

### Features

* frequencies: Include `--minimal-clade-size-to-estimate` command line option.
  [See PR 383](https://github.com/nextstrain/augur/pull/383)
* lbi: Include `--no-normalization` command line option.
  [See PR 380](https://github.com/nextstrain/augur/pull/380)

### Compatibility fixes

* export: Include `v1` subcommand to allow forwards compatibiliy with Augur v6 builds.
  [See PR 398](https://github.com/nextstrain/augur/pull/398)

### Bug fixes

* export: Include warning if using a mismatched v6 translate file.
  [See PR 392](https://github.com/nextstrain/augur/pull/392)
* frequencies: Fix determination of interval for clipping of non-informative pivots

## 5.3.0 (9 September 2019)

### Features

* export: Improve printing of error messages with missing or conflicting author
  data. [See issue 274](https://github.com/nextstrain/augur/issues/274)
* filter: Improve printing of dropped strains to include reasons why strains were
  dropped. [See PR 367](https://github.com/nextstrain/augur/pull/367)
* refine: Add support for command line flag `--keep-polytomies` to not resolve
  polytomies when producing a time tree.
  [See PR 345](https://github.com/nextstrain/augur/pull/345)

### Bug fixes

* Catch and throw error when there are duplicate strain names.
  [See PR 356](https://github.com/nextstrain/augur/pull/356)
* Fix missing annotation of "parent" attribute for the root node
* Run shell commands with more robust error checking.
  [See PR 350](https://github.com/nextstrain/augur/pull/350)
* Better handling of rerooting options for trees without temporal information.
  [See issue 348](https://github.com/nextstrain/augur/issues/348)

### Data

* Small fixes in geographic coordinate file

## 5.2.1 (4 August 2019)

### Bug fixes

* Print more useful error message if Python recursion limit is reached.
  [See issue 328](https://github.com/nextstrain/augur/issues/328)
* Print more useful error message if vcftools if missing.
  [See PR 312](https://github.com/nextstrain/augur/pull/321)

### Development

* Significantly relax version requirements specified in setup.py for biopython,
  pandas, etc... Additionally, move lesser used packages (cvxopt, matplotlib,
  seaborn) into an "extras_require" field. This should reduce conflicts with
  other pip installed packages.
  [See PR 323](https://github.com/nextstrain/augur/pull/323)

### Data

* Include additional country lat/longs in base data

## 5.2.0 (23 July 2019)

### Features

* ancestral: Adds a new flag `--output-sequences` and logic to support saving
  ancestral sequences and leaves from the given tree to a FASTA file. Also adds a
  redundant, more specific flag `--output-node-data` that will replace the current
  `--output` flag in the next major version release of augur. For now, we issue a
  deprecation warning when the `--output` flag is used. Note that FASTA output is
  only allowed for FASTA inputs and not for VCFs. We don't allow FASTA output for
  VCFs anywhere else and, if we did here, the output files would be very large.
  [See PR 293](https://github.com/nextstrain/augur/pull/293)

* frequencies: Allow `--method kde` flag to compute frequencies via KDE kernels.
  This complements existing method of `--method diffusion`. Generally, KDE
  frequencies should be more robust and faster to run, but will not project as
  well when forecasting frequencies into the future.
  [See PR 271](https://github.com/nextstrain/augur/pull/271)

### Bug fixes

* ancestral, traits, translate: Print warning if supplied tree is missing internal
  node names (normally provided by running `augur refine`).
  [See PR 283](https://github.com/nextstrain/augur/pull/283)

* Include pip in Conda enviroment file.
  [See PR 309](https://github.com/nextstrain/augur/pull/309)

### Documentation

* Document environment variables respected by Augur

### Development

* Remove matplotlib and seaborn from `setup.py` install. These are still called a
  few places in augur (like `titers.validate()`), but it was deemed rare enough
  that remove this from `setup.py` would ease general install for most users.
  Additionally, the ipdb debugger has been moved to dev dependencies.
  [See PR 291](https://github.com/nextstrain/augur/pull/291)

* Refactor logic to read trees from multiple formats into a function. Adds a new
  function `read_tree` to the `utils` module that tries to safely handle reading
  trees in multiple input formats.
  [See PR 310](https://github.com/nextstrain/augur/pull/310)


## 5.1.1 (1 July 2019)

### Features

* tree: Add support for the GTR+R10 substitution model.
* tree: Support parentheses in node names when using IQ-TREE.

### Bug fixes

* Use the center of the UK for its coordinates instead of London.
* filter: Mark `--output` required, which it always was but wasn't marked.
* filter: Avoid error when no excluded strains file is provided.
* export: Fix for preliminary version 2 schema support.
* refine: Correct error handling when the tree file is missing or empty.

### Documentation

* Add examples of Augur usage in the wild.
* Rename and reorganize CLI and Python API pages a little bit to make "where do
  I start learning to use Augur?" clearer to non-devs.

### Development

* Relax version requirements of pandas and seaborn.  The hope is this will make
  installation smoother (particularly alongside other packages which require
  newer pandas versions) while not encountering breaking changes in newer
  versions ourselves.


## 5.1.0 (29 May 2019)

### Documentation

* Documentation is now available online for the augur CLI and Python API via
  Read The Docs: <https://nextstrain-augur.readthedocs.io>.  The _latest_
  version on RTD points to the git master branch, and the _stable_ version to
  the most recent tagged release.  Instructions for building the docs locally
  are in the README.


## 5.0.0 (26 May 2019)

### Features

* ancestral: New option to `--keep-ambiguous`, which will not infer nucleotides at
  ambiguous (N) sites on tip sequences and instead leave as 'N'
  [See PR 280.](https://github.com/nextstrain/augur/pull/280)
* ancestral: New option to `--keep-overhangs`, which will not infer nucleotides for
  gaps on either side of the alignment and instead leave as '-'.
  [See PR 286.](https://github.com/nextstrain/augur/pull/286)
* clades: This module has been reconfigured to identify clade defining mutations on
  top of a reference rather than identifying mutations along the tree. The command
  line arguments are the same except for the addition of `--reference`, which
  explicitly passes in a reference sequence. If `--reference` is not defined, then
  reference will be drawn from the root node of the phylogeny by looking for
  `sequence` attribute attached to root node of `--tree`.
  [See PR 288.](https://github.com/nextstrain/augur/pull/288)
* refine: Revise rooting behavior. Previously `--root` took 'best', 'residual', 'rsq'
  and 'min_dev' as options. In this update `--root` takes 'best', least-squares',
  'min_dev' and 'oldest' as rooting options. This eliminates 'residual' and 'rsq'
  as options. This is a **backwards-incompatible** change. This requires updating
  TreeTime to version 0.5.4 or above.
  [See PR 263.](https://github.com/nextstrain/augur/pull/263)
* refine: Add `--keep-root` option that overrides `--root` specification to preserve
  tree rooting.
  [See PR 263.](https://github.com/nextstrain/augur/pull/263)
* refine: Add `--covariance` and `--no-covariance` options that specify TreeTime
  behavior.
  [See PR 263.](https://github.com/nextstrain/augur/pull/263)
* titers: This command now throws an `InsufficientDataException` if there are not
  sufficient titers to infer a model. This is paired with a new `--allow-empty-model`
  flag that proceeds past the `InsufficientDataException` and writes out a model
  JSON corresponding to an 'empty' model.
  [See PR 281.](https://github.com/nextstrain/augur/pull/281)
* By default JSONs are written with `index=1` to give a pretty-printed JSON. However,
  this adds significant file size to large tree JSONs. If the environment variable
  `AUGUR_MINIFY_JSON` is set then minified JSONs are printed instead. This mirror the
  explicit `--minify-json` argument available to `augur export`.
  [See PR 278.](https://github.com/nextstrain/augur/pull/278)

### Bug fixes

* export: Cast numeric values to strings for export.
  [See issue 287.](https://github.com/nextstrain/augur/issues/287)
* export: Legend order preserves ordering passed in by user for traits that have
  default colorings ('country' and 'region').
  [See PR 284.](https://github.com/nextstrain/augur/pull/284)
* refine: Previously, the `--root` argument was silently ignored when no timetree
  was inferred. Re-rooting with an outgroup is sensible even without a timetree.
  [See PR 282.](https://github.com/nextstrain/augur/pull/282)

## 4.0.0 (24 April 2019)

### Features

* distance: New interface for specifying distances between sequences. This is
  a **backwards-incompatible** change. Refer to `augur distance --help` for
  all the details.

* export: Add a `--minify-json` flag to omit indentation in Auspice JSONs.

### Bug fixes

* frequencies: Emit one-based coordinates (instead of zero-based) for KDE-based
  mutation frequencies

### Data

* Include additional country lat/longs in base data


## 3.1.8 (13 February 2019)

### Bug fixes

* titers: fix calculation of `mean_potentency` for model export

## 3.1.7 (5 February 2019)

### Bug fixes

* Update to TreeTime 0.5.3
* tree: Fix bug in printing causing errors in Python versions <3.6
* tree: Alter site masking to not be so memory intensive

## 3.1.6 (29 January 2019)

### Features

* filter: Allow negative matches to `--exclude-where`. For example,
 `--exclude-where country!=usa` would exclude all samples where metadata `country` does
 not equal `usa`.
* tree: Allow `--exclude-sites` to work with FASTA input. Ensure that indexing of input
 sites is one-based.

### Bug fixes

* fix loading of strains when loading titers from file, previously strains had not been
 filtered to match the tree appropriately

## 3.1.5 (13 January 2019)

### Features

* frequencies: Add `--ignore-char` and `--minimal-clade-size` as options.
* frequencies: Include `--stiffness` and `--inertia` as options.
* titers: Allow multiple titer date files in `--titers` import.

### Bug fixes

* filter: Fix `--non-nucleotide` call to include `?` as allowed character.
* tree: Fix `--method raxml` to properly delimit interim RAxML output so that
  simultaneous builds don't conflict.

### Data

* Include additional country lat/longs in base data

## 3.1.4 (1 January 2019)

### Bug fixes

* frequencies: Include `counts` in `augur frequencies` output JSON to support
  downstream plotting.

### Data

* Include additional country lat/longs in base data

## 3.1.3 (29 December 2018)

### Features

* filter: Add `--non-nucleotide` option to remove sequences with non-conforming
  nucleotide characters.

### Bug fixes

* Revise treatment of `-`, ` ` in `augur parse` to leave `-` as is and remove white
  space. Also delimit `[` and `]` to `_`.
* Fix bug in naming of temp IQTREE fixes to prevent conflicts from simultaneous builds.

### Data

* Include additional country lat/longs in base data

### Development

* Remove non-modular measles build in favor of [nextstrain/measles](https://github.com/nextstrain/measles) repo.

## 3.1.2 (21 December 2018)

### Bug fixes

* Update dependencies

## 3.1.1 (21 December 2018)

### Bug fixes

* filter: Fix `--include-where`. Adds an `all_seq` variable needed by the logic to
  include records by value. This was previously working for VCF but threw an exception
  for sequences in FASTA format.
* Update flu reference viruses and lat longs.
* Update dependencies

## 3.1.0 (18 December 2018)

### Features

* reconstruct-sequences: Include `augur reconstruct-sequences` module that reconstructs
  alignments from mutations inferred on the tree
* distance: Include `augur distance` module that calculates the distance between amino
  acid sequences across entire genes or at a predefined subset of sites
* lbi: Include `augur lbi` module that calculates local branching index (LBI) for a
  given tree and one or more sets of parameters.
* frequencies: Include `--method kde` as option to `augur frequencies`, separate from the
  existing `--method diffusion` logic. KDE frequencies are faster and better for smaller
  clades but don't extrapolate as well as diffusion frequencies.
* titers: Enable annotation of nodes in a tree from the substitution model

## 3.0.5.dev1 (26 November 2018)

### Bug fixes

* translate: Nucleotide ("nuc") annotation for non-bacterial builds starts at 0
  again, not 1, fixing a regression.

### Documentation

* Schemas: Correct coordinate system description for genome start/end
  annotations.


## 3.0.4.dev1 (26 November 2018)

### Bug fixes

* validate: Fix regression for gene names containing an asterisk.

### Development

* Fix Travis CI tests which were silently not running.


## 3.0.3.dev1 (26 November 2018)

### Features

* refine: Add a `--clock-std-dev` option

* traits: Add a `--sampling-bias-correction` option for mugration model

* validate: Gene names in tree annotations may now contain hyphens.  Compatible
  with Auspice version 1.33.0 and later.

* All JSON is now emitted with sorted keys, making it easier to diff and run
  other textual comparisons against output.

### Bug fixes

* filter: Only consider A, T, C, and G when calculating sequence length for the
  `--min-length` option.

* filter: Allow comments in files passed to `--exclude`.

* filter: Ignore case when matching trait values against excluded values.

* Normalize custom geographic names to lower case for consistent matching.

### Data

* Fix typo in geographic entry for `netherlands`.

* Schemas: Reconcile naming patterns used in gene definitions and tree
  annotations.

### Development

* Upgrade TreeTime dependency to 0.5.x and at least 0.5.1.

* Add an `environment.yml` file for use with `conda env create`.

* Stop testing under Python 2.7 on Travis CI.


## 3.0.2.dev1 (27 September 2018)

### Bug fixes

* translate: Fix broken `--help` message


## 3.0.1.dev1 (27 September 2018)

### Features

* align and tree: The --nthreads option now accepts the special value "auto" to
  automatically set the number of threads to the number of CPU cores available.

* Alias `augur --version` to `augur version`

### Bug fixes

* tree: The --nthreads option is now respected.  Previously all tree builders
  were ignoring the value and using either 2 threads (RAxML, IQ-TREE) or as
  many threads as cores (FastTree, if the OpenMP version).

* translate: Check for and, if necessary pad, nucleotide sequences which aren't
  a multiple of 3 earlier to avoid errors later.

* export: Optionally write inferred nucleotide and amino acid sequences (or
  mutations) to a separate file.

* export: Omit genes with no amino acid mutations.

* validate: Allow underscores in gene names.

* refine: Remove unused --nthreads argument.

* ancestral, filter, tree, refine: Exit 1 instead of -1 on error.

* Print the help message, instead of throwing an exception, when `augur` is run
  without arguments.

### Documentation

* Briefly describe each command in its `--help` output and in the global `augur
  --help` output.

* Revamp README to emphasize new, modular augur and make it suitable for
  inclusion on PyPi.

* Reconciled conflicting license declarations; augur is AGPLv3 (not MIT)
  licensed like the rest of Nextstrain.

* Include URLs for bug reports, the change log, and the source on PyPi.

### Data

* Geographic coordinates added for the Netherlands and the Philippines.

### Development

* Reset the `release` branch when rewinding a failed local release process.

* Refactor the augur program and command architecture for improved
  maintainability.


## 3.0.0.dev3 (4 September 2018)

### Development

* Use an allowed Topic classifier so we can upload to PyPi

* Ignore distribution egg-info build files


## 3.0.0.dev2 (4 September 2018)

### Features

* Export: Add safety checks for optional annotations and geo data

* Include more lat/longs in the default geo data

### Development

* Add release tooling

* Document the release process and a few development practices

* Travis CI: Switch to rebuilding the Docker image only for new releases

* Remove ebola, lassa, tb, WNV, and zika builds now in their own repos.  These
  builds are now available at URLs like <https://github.com/nextstrain/ebola>,
  for example.


## 3.0.0.dev1 (unreleased)

### Development

* Start versioning augur beginning with 3.0.0.  A new `augur version` command
  reports the running version.
