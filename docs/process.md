# Process


### how to run
The process scripts (and the Process class) is designed to analyse one single dataset (i.e. not multiple segments), coming from one single JSON.
Similar to _Prepare_, it is designed off a `config` dict that should be all that is needed to run the analysis.
For bespoke analysis, a new class may be created which inherits from `Process`.
Since the config files across segments may be identical but for the filenames, these can be added as command line arguments on a dataset-by-dataset basis (e.g. see flu).
Most stages will save their output which will be restored (if valid) for subsequent analysis.

### input
Normally, everything is done off a JSON via [prepare](./prepare.md).
Titer information is loaded here, although this is a work in progress.

### analysis components
* align
* [estimate mutation frequencies](./frequencies.md)
* [build tree (RAxML)](./phylogenies.md)
* [TimeTree](./phylogenies.md) - temporal filtering, ancestral reconstruction + node dating.
* geographical inference
* [estimate tree frequencies](./frequencies.md)
* Titer models
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
  * `controls` {dict} (default: `{}`)
    * `geographic_location`
    * `authors`
