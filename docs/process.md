# Process

The process scripts (and the Process class) are designed to analyse data one prepared JSON at a time.
This allows segmented viruses, or multiple lineages to be analysed independently.
Similar to [_prepare_](prepare.md), it is designed off a `config` dict that should be all that is needed to run the analysis.
For bespoke analysis, a new class may be created which inherits from `Process`.
Most scripts (see [zika](../zika/zika.process.py) and [flu](../flu/flu.process.py)) use command line arguments to dynamically change the `config` dictionary.
A number of intermediate files (e.g. trees, alignments) are created as well as JSONs to be visualized in Auspice.
While _process_ can be computationally expensive, as long as the underlying _prepared JSON_ is unchanged the script should be able to restore itself when you rerun it.

### analysis components
* align
* [estimate mutation frequencies](./frequencies.md)
* [build tree (RAxML)](./phylogenies.md)
* [TimeTree](./phylogenies.md) - temporal filtering, ancestral reconstruction + node dating.
* geographical inference
* [estimate tree frequencies](./frequencies.md)
* [annotate tree with scores](./scores.md) - calculate sequence- and metadata-based scores for each node in the tree
* Titer models (See [Neher et al, PNAS, 2016](http://www.pnas.org/content/113/12/E1701.abstract) )
* Identify predetermined clades
* [Export for auspice](./auspice_output.md)

### Config dict

##### Basic settings
* `dir` the current directory - not _augur_ but the pathogen itself
* `output` dict with 2 keys: `data` & `auspice`. (Default: `"output": {"data": "processed","auspice": "auspice",}`
  * `data`: progress files intermediate FASTAs, nexus trees etc
  * `auspice`: the finished JSONs for display in auspice
* `in` {string}: The prepared JSON file. Sometimes defined by command line arguments
* `clean` {bool} (default `False`) Run a clean build by removing any (previously created) files
* `subprocess_verbosity_level` {int} (default `0`). Control the amount of information displayed during alignment and tree building. Higher numbers -> more output.


##### Analysis settings
* `geo_inference` {`False` || array of strings} (Default: `False`) what traits to perform geographic inference (mugration model) upon
* `geo_inference_options`: {dict}. Keys:
  * `confidence` {bool} (default: `True`) Include (normalized) likelihoods from any geographic inference analysis.
  * `root_state` {dict} (see zika)
* `temporal_confidence` {bool} (Default: `True`) Include 90% (normalized marginal likelihood) confidence intervals from TreeTime dating analysis.
* `estimate_mutation_frequencies` {bool} (default: `False`)
* `pivot_spacing` {float} (default: not present) necessary if you want to calculate pivots through `get_pivots_via_spacing`. See [Frequencies](./frequencies.md)
* `newick_tree_options` {dict} (default: `{}`) See [Phylogenies](./phylogenies.md)
* `clock_filter` {`False` | dict} See [Phylogenies](./phylogenies.md)
* `timetree_options` {dict} See [Phylogenies](./phylogenies.md)
* `epitope_mask` {string} a tab-delimited file containing epitope mask names in the first column and a bitmask of sites in the protein coding sequence of a virus segment associated with epitopes
* `epitope_mask_version` {string} the name of an epitope mask defined in the `epitope_mask` file to use for analyses
* `predictors` {array of strings or dict} a list of attributes annotated to each tree node to use for the fitness model (e.g., `"['ep']"` or `"['cTiter', 'ep']"`) or a dictionary of predictors and their corresponding precalculated model parameters and global standard deviations (e.g., `{'ep': [0.33, 1.31]}`)

##### Auspice output settings
* `auspice`
  * `panels` {array} (default: `['tree', 'map', 'entropy']`)
  * `extra_attr` {array} (default: `[]`)
  * `date_range`
  * `analysisSlider`: optional. If specified, this should be a key present in `color_options`, with `type`: continuous
  * `color_options` (default: `"num_date": {...}, "gt": {...}`)
    * `<trait name>`
      * `key`
      *  `legendTitle`: String appearing in Auspice legend
      * `menuItem`: String appearing in Auspice colorBy drop down menu
      * `type`: continuous / integer / discrete
  * `controls` {dict} (default: `{}`) **THIS IS DEPRECATED AND WILL BE REMOVED SOON**
    * `geographic_location`
    * `authors`
  * `defaults`: {dict} Auspice defaults (exported as-is to the `meta.json` file)
    * Keys (and the types of their values) may include `colorBy` (string), `geoResolution` (string), `distanceMeasure` (string), `mapTriplicate` (bool)
