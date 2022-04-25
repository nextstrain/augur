# CHANGELOG

## __NEXT__


## 15.0.0 (15 April 2022)

### Major Changes

* export: Move extensions block to meta [#888][] (@corneliusroemer)

    Note: this is technically a breaking change, but the misplaced extensions block was added in version 14.1.0 and intended for internal use by Nextclade. We don't expect any users to be impacted by this.

### Features

* filter: Support relative dates for `--min-date` and `--max-date` [#740][] (@benjaminotter, @victorlin)
* frequencies: Support relative dates for `--min-date` and `--max-date` [#890][] (@victorlin, @huddlej)
* docs: Add page describing how to use your own tree builder [#884][] (@jameshadfield, @joverlee521)

### Bug Fixes

* parse: Fix `--fix-dates` [#878][] (@huddlej)
* tree: Fix internal comment on excluded sites [#880][] (@jameshadfield)

[#888]: https://github.com/nextstrain/augur/pull/888
[#740]: https://github.com/nextstrain/augur/pull/740
[#890]: https://github.com/nextstrain/augur/pull/890
[#884]: https://github.com/nextstrain/augur/pull/884
[#878]: https://github.com/nextstrain/augur/pull/878
[#880]: https://github.com/nextstrain/augur/pull/880

## 14.1.0 (31 March 2022)

### Features

* schemas: Extend export v2 schema to support an array of trees [#851][] (@tsibley)
* schemas: Add JSON schemas for our root-sequence and tip-frequencies sidecars [#852][] (@tsibley)
* schemas: Add JSON schema for measurements sidecar [#859][] (@joverlee521)
* filter: Send warnings to stderr to be consistent with other warnings [#862][] (@victorlin)
* export: Allow an extensions block in auspice config & dataset JSONs [#865][] (@jameshadfield)
* export: Allow skipping of input/output schema validation [#865][] (@jameshadfield)
* export:  Order keys in dataset for easier reading [#868][] (@jameshadfield)

### Bug Fixes

* parse: Fix typo in internal variable name [#850][] (@emmahodcroft)

[#850]: https://github.com/nextstrain/augur/pull/850
[#851]: https://github.com/nextstrain/augur/pull/851
[#852]: https://github.com/nextstrain/augur/pull/852
[#859]: https://github.com/nextstrain/augur/pull/859
[#862]: https://github.com/nextstrain/augur/pull/862
[#865]: https://github.com/nextstrain/augur/pull/865
[#868]: https://github.com/nextstrain/augur/pull/868

## 14.0.0 (8 February 2022)

### Major Changes

* Drop support for Python 3.6, add support for 3.9 and 3.10 [#822][] (@victorlin)

### Features

* refine: Enable bootstrap support by passing confidence values through to Auspice JSONs from `augur refine` node data JSONs [#839][] (@huddlej)
* tree: Allow users to override default tree builder arguments with a new `--override-default-args` flag [#839][] (@huddlej)
* clades: Allow descendant clades to be defined by explicitly inheriting from ancestral clade names [#846][] (@corneliusroemer)

### Bug Fixes

* tree: Fix segmentation fault that can occur when user-provided tree builder args conflict hardcoded defaults for IQ-TREE's. The new `--override-default-args` flag allows users to override the defaults that conflict with their values. [#839][] (@huddlej)
* filter/utils: fix year-only and numeric date handling [#841][] (@victorlin)
* CI: test earliest supported Biopython versions in matrix, remove redundant installs [#843][] (@victorlin)

[#822]: https://github.com/nextstrain/augur/pull/822
[#841]: https://github.com/nextstrain/augur/pull/841
[#843]: https://github.com/nextstrain/augur/pull/843
[#846]: https://github.com/nextstrain/augur/pull/846
[#839]: https://github.com/nextstrain/augur/pull/839

## 13.1.2 (28 January 2022)

### Features

### Bug Fixes

* filter: Set default filter priority to negative infinity instead of 0. [#835][] (@huddlej)
* filter: Cast strains to list when subsetting metadata to prevent pandas FutureWarning. [#833][] (@victorlin)
* export v2: Allow periods (.) in accession numbers. [#832][] (@tsibley)

[#835]: https://github.com/nextstrain/augur/pull/835
[#833]: https://github.com/nextstrain/augur/pull/833
[#832]: https://github.com/nextstrain/augur/pull/832

## 13.1.1 (21 January 2022)

### Bug Fixes

* dependencies: Relax upper bounds on development dependencies. [#821][] and [#830][] (@victorlin, @huddlej)
* refine: Inform users when their requested root sequence is not present in the alignment. [#816][] (@corneliusroemer)
* utils: Fix crash with `read_metadata` when numexpr is installed. [#826][] (@victorlin)

[#816]: https://github.com/nextstrain/augur/pull/816
[#821]: https://github.com/nextstrain/augur/pull/821
[#826]: https://github.com/nextstrain/augur/pull/826
[#830]: https://github.com/nextstrain/augur/pull/830

## 13.1.0 (10 December 2021)

### Features

* schemas: Add "$id" key to Auspice config schemas so we have a way of referring to these. [#806][] (@tsibley)

### Bug Fixes

* filter: Fix groupby with incomplete dates. [#808][] (@victorlin)

[#806]: https://github.com/nextstrain/augur/pull/806
[#808]: https://github.com/nextstrain/augur/pull/808

## 13.0.4 (8 December 2021)

### Bug Fixes

* dependencies: Replace deprecated mutable sequence interface for BioPython. [#788][] (@Carlosbogo)
* dependencies: Fix backward compatibility with BioPython. [#801][] (@huddlej)
* data: Add latitude and longitude details for "Reunion". [#791][] (@corneliusroemer)
* filter: Use pandas functions to determine subsample groups. [#794][] and [#797][] (@victorlin)
* filter: Add clarity to help message and output of probabilistic sampling. [#792][] (@victorlin)

[#788]: https://github.com/nextstrain/augur/pull/788
[#791]: https://github.com/nextstrain/augur/pull/791
[#792]: https://github.com/nextstrain/augur/pull/792
[#794]: https://github.com/nextstrain/augur/pull/794
[#797]: https://github.com/nextstrain/augur/pull/797
[#801]: https://github.com/nextstrain/augur/pull/801

## 13.0.3 (19 November 2021)

### Bug Fixes

* tree: Handle compressed alignment when excluding sites. [#786][] (@huddlej)
* docs: Fix typos ([ce0834c][]) and clarify exclude sites inputs ([5ad1574][]). (@corneliusroemer)

[#786]: https://github.com/nextstrain/augur/pull/786
[ce0834c]: https://github.com/nextstrain/augur/commit/ce0834c476abc9ee99785fa930608218b7d78990
[5ad1574]: https://github.com/nextstrain/augur/commit/5ad157485015623883c6b637d247459f906b63cb

## 13.0.2 (12 October 2021)

### Bug Fixes

* dependencies: Support latest versions of BioPython. [#777][] (@huddlej)
* tree: Allow users to specify arbitrary IQ-TREE models. [#776][] (@huddlej)

[#776]: https://github.com/nextstrain/augur/pull/776
[#777]: https://github.com/nextstrain/augur/pull/777

## 13.0.1 (1 October 2021)

### Bug Fixes

* docs: Fix broken link to latitude/longitude documentation. [#766][] (@victorlin)
* filter: Fix reproducibility of subsampling by using the user-defined random seed in all random function calls and by ordering strain sets as lists prior to adding strains to group-by priority queues. [#772][] (@huddlej)

[#766]: https://github.com/nextstrain/augur/pull/766
[#772]: https://github.com/nextstrain/augur/pull/772

## 13.0.0 (17 August 2021)

### Major Changes

* filter: Skip metadata records with ambiguous month information in the `date` column when grouping by month instead of randomly generating month values for those records. This change alters the behavior of the `filter` command for metadata with ambiguous month values. For these data, consider using `--group-by year` instead of `--group-by year month`. [#761][] (@huddlej)

### Features

* filter: When grouping by year or month, report the number of strains skipped due to ambiguous year and month both in the summary report at the end of filtering and in the `--output-log` contents [#761][] (@huddlej)

[#761]: https://github.com/nextstrain/augur/pull/761

## 12.1.1 (13 August 2021)

### Bug Fixes

* filter: Fix parsing of missing data in metadata [#758][] (@huddlej)
* filter: Fix probabilistic sampling with small values [#759][] (@huddlej)

[#758]: https://github.com/nextstrain/augur/pull/758
[#759]: https://github.com/nextstrain/augur/pull/759

## 12.1.0 (12 August 2021)

### Features

* export: Add support for custom legend and color scale specifications in Auspice config files [#727][] (@jameshadfield)
* utils: Add support for compressed strain name files (e.g., "include.txt.gz") [#730][] (@benjaminotter)
* filter: Rewrite internal logic to use pandas DataFrames ([#743][]), define filters and subsampling logic as individual functions ([#745][] and [#746][]), and iterate through chunks of metadata instead of loading all records into memory at once ([#750][]) (@tsibley, @huddlej)

### Bug Fixes

* distance: Change numeric type of distance output to float [#729][] (@benjaminotter)
* filter: Disable probabilistic sampling when users provide `--sequences-per-group` [#737][] (@benjaminotter)
* export: Provide correct missing file error messages for metadata and node data JSON inputs [#752][] (@benjaminotter)

[#727]: https://github.com/nextstrain/augur/pull/727
[#729]: https://github.com/nextstrain/augur/pull/729
[#730]: https://github.com/nextstrain/augur/pull/730
[#737]: https://github.com/nextstrain/augur/pull/737
[#743]: https://github.com/nextstrain/augur/pull/743
[#745]: https://github.com/nextstrain/augur/pull/745
[#746]: https://github.com/nextstrain/augur/pull/746
[#750]: https://github.com/nextstrain/augur/pull/750
[#752]: https://github.com/nextstrain/augur/pull/752

## 12.0.0 (13 April 2021)

### Major Changes

* filter: Date bounds (`--min-date` and `--max-date`) are now inclusive instead of exclusive such that records matching the given dates will pass date filters [#708][] (@benjaminotter)

### Bug Fixes

* refine: Recommend an alternate action when skyline optimization fails [#712][] (@huddlej)

### Features

* distance: Count insertion/deletion events once in pairwise distances [#698][] (@huddlej, @benjaminotter)
* distance: Optionally ignore specific list of characters defined in a distance map's top-level `ignored_characters` list [#707][] (@benjaminotter)
* filter: Allow `--subsample-max-sequences` without `--group-by` [#710][] (@benjaminotter)
* tree: Prefer `iqtree2` binary over `iqtree` when possible [#711][] (@benjaminotter)

[#698]: https://github.com/nextstrain/augur/pull/698
[#707]: https://github.com/nextstrain/augur/pull/707
[#708]: https://github.com/nextstrain/augur/pull/708
[#710]: https://github.com/nextstrain/augur/pull/710
[#711]: https://github.com/nextstrain/augur/pull/711
[#712]: https://github.com/nextstrain/augur/pull/712

## 11.3.0 (19 March 2021)

### Bug Fixes

* filter: Clarify how the `--priority` input affects subsampling in the command line help documentation [#695][]
* tests: Clean up outputs created by tests [#703][], ignore log files [#701][], and remove unnecessary Conda environment file [#702][]

### Features

* io: Add new `io` module with `open_file`, `read_sequences`, and `write_sequences` functions that support compressed inputs and outputs [#652][]
* parse, index, filter, mask: Add support for compressed inputs/outputs [#652][]
* export v2: Add optional `data_provenance` field to auspice JSON output for better provenance reporting in Auspice [#705][]

[#652]: https://github.com/nextstrain/augur/pull/652
[#695]: https://github.com/nextstrain/augur/pull/695
[#701]: https://github.com/nextstrain/augur/pull/701
[#702]: https://github.com/nextstrain/augur/pull/702
[#703]: https://github.com/nextstrain/augur/pull/703
[#705]: https://github.com/nextstrain/augur/pull/705

## 11.2.0 (8 March 2021)

### Bug Fixes

* ancestral: Mask positions that are ambiguous in all tip sequences before inferring ancestral sequence states, to avoid assigning arbitrary ancestral values based on rounding errors [#682][]
* titers: Add missing `kwargs` attribute to `TiterCollection` class [#690][]

### Documentation

* Update API documentation to include newer Python modules and the `index` subcommand [#687][]
* Remove Zika and TB tutorials in favor of copies in docs.nextstrain.org [#689][]

### Features

* filter: Enable filtering by metadata only such that sequence inputs/outputs are optional and metadata/strain list outputs are now possible [#679][]
* filter: Enable extraction of sequences from multiple lists of strains with a new `--exclude-all` flag and support for multiple inputs to the `--include` argument [#679][]

[#679]: https://github.com/nextstrain/augur/pull/679
[#682]: https://github.com/nextstrain/augur/pull/682
[#687]: https://github.com/nextstrain/augur/pull/687
[#689]: https://github.com/nextstrain/augur/pull/689
[#690]: https://github.com/nextstrain/augur/pull/690

## 11.1.2 (16 February 2021)

### Bug Fixes

* index: Remove call to deprecated BioPython SeqIO.close method [#684][]

[#684]: https://github.com/nextstrain/augur/pull/684

## 11.1.1 (16 February 2021)

### Bug Fixes

* filter: Retry probabilistic subsampling when it doesn't select any samples [#676][]
* titers: Skip infinite log titer values caused by zero-valued raw titers [#677][]
* setup: Include license file with distribution artifacts instead of Python installation root [#678][]

[#676]: https://github.com/nextstrain/augur/pull/676
[#677]: https://github.com/nextstrain/augur/pull/677
[#678]: https://github.com/nextstrain/augur/pull/678

## 11.1.0 (12 February 2021)

### Bug Fixes

* align/utils: Improve external command errors [#638][]
* filter: Fix parsing of priorities files to allow spaces in sequence IDs [#668][]
* utils: Ensure columns `strain` and `name` in metadata get parsed as strings [#669][]

### Features

* index/filter: Add new `index` subcommand and optional `--sequence-index` argument for `filter` command to enable filtering without inspecting sequences [#651][]
* titers: Bump supported cvxopt version to latest 1.* release [#671][]

[#638]: https://github.com/nextstrain/augur/pull/638
[#651]: https://github.com/nextstrain/augur/pull/651
[#668]: https://github.com/nextstrain/augur/pull/668
[#669]: https://github.com/nextstrain/augur/pull/669
[#671]: https://github.com/nextstrain/augur/pull/671

## 11.0.0 (22 January 2021)

### Major Changes

* filter: Use probabilistic sampling by default when requesting a maximum number of sequences to subsample with `--subsample-max-sequences`. Adds `--no-probabilistic-sampling` flag to disable this default behavior and prevent users from requesting fewer maximum sequences than there are subsampling groups. [#659][]

[#659]: https://github.com/nextstrain/augur/pull/659

## 10.3.0 (14 January 2021)

### Bug Fixes

* scripts: Fix typo in `verify_meta_json.py` [#656][] (@felixonmars)
* CI: Use the expected Python version in conda environments [#658][]
* CI: Minimize codecov feedback [#661][]

### Features

* frequencies: Add `--pivot-interval-units` argument with support for weekly pivots [#660][]
* frequencies: Add support for ISO dates for `--min-date` and `--max-date` arguments [#660][]

[#656]: https://github.com/nextstrain/augur/pull/656
[#658]: https://github.com/nextstrain/augur/pull/658
[#660]: https://github.com/nextstrain/augur/pull/660
[#661]: https://github.com/nextstrain/augur/pull/661

## 10.2.0 (1 January 2021)

### Features

* filter: Add `--probablistic-sampling` flag to allow subsampling with `--subsample-max-sequences` when the number of groups exceeds the requested number of samples [#629][]
* scripts: Add script to identify emerging clades from existing Nextstrain build JSONs [#653][]
* docs: Add instructions to update conda installations prior to installing Augur [#655][]

[#629]: https://github.com/nextstrain/augur/pull/629
[#653]: https://github.com/nextstrain/augur/pull/653
[#655]: https://github.com/nextstrain/augur/pull/655

## 10.1.1 (16 November 2020)

### Bug Fixes

* dependencies: Require the most recent minor versions of TreeTime (0.8.X) to fix numpy matrix errors [#633][]

[#633]: https://github.com/nextstrain/augur/pull/633

## 10.1.0 (13 November 2020)

### Features

* docs: Migrate non-reference documentation to docs.nextstrain.org [#620][]
* filter: Add `--exclude-ambiguous-dates-by` flag to enable exclusion of samples with ambiguous dates [#623][] and [#631][]

[#620]: https://github.com/nextstrain/augur/pull/620
[#623]: https://github.com/nextstrain/augur/pull/623
[#631]: https://github.com/nextstrain/augur/pull/631

## 10.0.4 (6 November 2020)

### Bug Fixes

* tree: Use a more generic approach to escape special characters from alignment sequence names prior to running IQ-TREE [#625][]
* filter: Reduce memory usage by not reading sequences into memory [#627][]

[#625]: https://github.com/nextstrain/augur/pull/625
[#627]: https://github.com/nextstrain/augur/pull/627

## 10.0.3 (23 October 2020)

### Bug Fixes

* refine: Report divergence by number of mutations as an integer instead of a floating point value [#618][]
* validate: Allow internal nodes with a single child and do not allow duplicate node names [#621][]

[#618]: https://github.com/nextstrain/augur/pull/618
[#621]: https://github.com/nextstrain/augur/pull/621

## 10.0.2 (8 September 2020)

### Bug Fixes

* align: Remove references to BioPython's deprecated `Alphabet` attributes [#615][]
* Pin BioPython dependency to a max supported version to prevent breaking changes to augur in the future [#615][]

[#615]: https://github.com/nextstrain/augur/pull/615

## 10.0.1 (8 September 2020)

### Bug Fixes

* ancestral: Clarify default values for inference of ambiguous bases [#613][]

[#613]: https://github.com/nextstrain/augur/pull/613

## 10.0.0 (17 August 2020)

### Major Changes

* Remove Snakemake as a dependency of the augur Python package [#557][]
  * Updated documentation to reflect external Snakemake dependency [#600][] and Snakemake's required `--cores` argument [#599][]
* utils: `read_colors` refactor [#588][]
  * raises an exception when the requested color file is missing instead of printing a warning to stdout
  * splits out logic to parse colors file into separate classes (`util_support/color_parser.py` and `util_support/color_parser_line.py`) with unit tests
* utils: `read_metadata` interface improvements
  * raises exceptions when 1) input file is missing or unreadable or 2) required columns (`strain` or `name`) are missing instead of failing silently [#584][]
  * automatically detects delimiter in metadata file instead of assuming delimiter based on filename extension [#587][]
* utils: `read_node_data` interface improvements [#595][], [#605][]
  * exits with a nonzero code when node data node names don't match tree nodes and when the input tree cannot be loaded
  * refactors logic to read node data into separate classes with unit tests

### Bug Fixes

* ancestral: Fix docstring for `collect_mutations_and_sequences` [4c474a9][]
* parse: Fix date parsing bug caused by a change in the API for `parse_time_string` in pandas 1.1.0 [#601][]
* refine: Enable divergence unit scaling without timetree [e9b3eec][]
* tree: Use IQ-TREE's `-nt AUTO` mode when users request more threads than there are input sequences, avoiding an IQ-TREE error [#598][]

### Features

* docs: Document support for installation from Bioconda [#604][]
* filter: Add `--subsample-max-sequences` argument to limit the maximum number of sequences to be included in subsampled output [#593][]
* mask: Add `--mask-invalid` flag to mask invalid nucleotides from FASTA files [#592][]

[#557]: https://github.com/nextstrain/augur/pull/557
[#584]: https://github.com/nextstrain/augur/pull/584
[#587]: https://github.com/nextstrain/augur/pull/587
[#588]: https://github.com/nextstrain/augur/pull/588
[#592]: https://github.com/nextstrain/augur/pull/592
[#593]: https://github.com/nextstrain/augur/pull/593
[#595]: https://github.com/nextstrain/augur/pull/595
[#598]: https://github.com/nextstrain/augur/pull/598
[#599]: https://github.com/nextstrain/augur/pull/599
[#600]: https://github.com/nextstrain/augur/pull/600
[#601]: https://github.com/nextstrain/augur/pull/601
[#604]: https://github.com/nextstrain/augur/pull/604
[#605]: https://github.com/nextstrain/augur/pull/605
[e9b3eec]: https://github.com/nextstrain/augur/commit/e9b3eec670b9603874e195cc1ccd4f3c1aeef5dd
[4c474a9]: https://github.com/nextstrain/augur/commit/4c474a96232e9cc333e3fc4c0971336a090b703c

## 9.0.0 (29 June 2020)

### Major Changes

* align: The API to the `read_sequences` function now returns a list of sequences instead of a dictionary [#536][]

### Bug Fixes

* align: Prevent duplicate strains warning when using `--reference-name` [#536][]
* docs: Sync and deduplicate installation documentation from README to main docs [#578][]
* export: Flexibly disambiguate multiple publications by the same author [#581][]
* frequencies: Avoid interpolation of a single data point during frequency estimation with sparse data [#569][]
* parse: Actually remove commas during prettify when this behavior is requested [#573][]
* tests: Always use the local helper script (`bin/augur`) to run tests instead of any globally installed augur executables [#527][]
* tree: Keep log files after trees are built [#572][]
* utils: Do not attempt to parse dates with only ambiguous months (e.g., 2020-XX-01) [#532][]
* utils: Parse `name` column of metadata as a data field instead of a pandas DataFrame attribute [#564][]

### Features

* docs: Updates description of how missing data are handled by `augur traits`
* filter: Add support for ISO 8601 dates (YYYY-MM-DD) for `--min-date` and `--max-date` [#568][]
* tests: Add tests for utilities (ambiguous date parsing [#532][] and `run_shell_command` [#577][]), parse [#573][], and translate [#546][]
* tree: Allow VCF input without an `--exclude-sites` argument [#565][]

[#527]: https://github.com/nextstrain/augur/pull/527
[#532]: https://github.com/nextstrain/augur/pull/532
[#536]: https://github.com/nextstrain/augur/pull/536
[#546]: https://github.com/nextstrain/augur/pull/546
[#564]: https://github.com/nextstrain/augur/pull/564
[#565]: https://github.com/nextstrain/augur/pull/565
[#568]: https://github.com/nextstrain/augur/pull/568
[#569]: https://github.com/nextstrain/augur/pull/569
[#572]: https://github.com/nextstrain/augur/pull/572
[#573]: https://github.com/nextstrain/augur/pull/573
[#577]: https://github.com/nextstrain/augur/pull/577
[#578]: https://github.com/nextstrain/augur/pull/578
[#581]: https://github.com/nextstrain/augur/pull/581

## 8.0.0 (8 June 2020)

### Major Changes

* utils: Add a consolidated generic `load_mask_sites` function and specific `read_mask_file` and `read_bed_file` functions for reading masking sites from files. Changes the Python API by moving mask-loading functionality out of augur mask and tree into utils [#514][] and [#550][]
* mask: Parse BED files as zero-indexed, half-open intervals [#512][]

### Bug Fixes

* traits: Export mugration models to the same output directory as traits JSON [#544][]
* Explicitly open files with UTF-8 file encoding [#499][], [#503][], and [#560][]
* refine: Only request confidence intervals from TreeTime when no clock rate is provided [#548][]
* refine: Catch failed skyline optimization [#556][]

### Features

* align: Report insertions stripped during alignment [#449][]
* Require minimum pandas version of 1.0.0 [#488][]
* parse: Reduce memory use and clarify code with standard Python idioms [#496][]
* mask: Allow masking of specific sites passed by the user with `--mask-sites` and masking of a fixed number of sites from the beginning or end of each sequence with `--mask-from-beginning` and `--mask-from-end` [#512][]
* clades, import: Use `defaultdict` to simplify code [#533][]
* tests: Add initial functional tests of the augur command line interface using Cram [#542][]
* refine: Add a `--seed` argument to set the random seed for more reproducible outputs across runs [#542][]
* ancestral, refine, and traits: Print the version of TreeTime being used for these commands [#552][]
* filter: Add support for flexible pandas-style queries with new `--query` argument [#555][]
* export: Allow display defaults for transmission lines [#561][]

[#449]: https://github.com/nextstrain/augur/pull/449
[#488]: https://github.com/nextstrain/augur/pull/488
[#496]: https://github.com/nextstrain/augur/pull/496
[#499]: https://github.com/nextstrain/augur/pull/499
[#503]: https://github.com/nextstrain/augur/pull/503
[#512]: https://github.com/nextstrain/augur/pull/512
[#514]: https://github.com/nextstrain/augur/pull/514
[#533]: https://github.com/nextstrain/augur/pull/533
[#542]: https://github.com/nextstrain/augur/pull/542
[#544]: https://github.com/nextstrain/augur/pull/544
[#548]: https://github.com/nextstrain/augur/pull/548
[#550]: https://github.com/nextstrain/augur/pull/550
[#552]: https://github.com/nextstrain/augur/pull/552
[#555]: https://github.com/nextstrain/augur/pull/555
[#556]: https://github.com/nextstrain/augur/pull/556
[#560]: https://github.com/nextstrain/augur/pull/560
[#561]: https://github.com/nextstrain/augur/pull/561

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
