# Prepare

Prepare takes fauna data as input and applies filtering, subsampling and checking of metadata to produce a single JSON (per segment) which can be analysed by [_process_](process.md).
Prepare scripts define a config dictionary which typically contains all the information needed.
If you need something a bit more bespoke that can't be done through the config, the idea would be to create a new class which inherits from `Prepare` and use that instead.

For examples see: [zika](../zika/zika.prepare.md) (simple) and [flu](../flu/flu.prepare.md) (more complicated).

### `Config` dictionary
The `config` dict is passed to the `Prepare` class.
Many options have defaults provided, and if so are not required to be in the config file.

#### general settings (including some auspice-specific settings)
* `dir`: the current directory - not _augur_ but the virus itself
* `file_prefix`: string used to name the JSONs
* `title`: string used for display in auspice (optional - `file_prefix` used if not specified)
* `maintainer`: array containing two strings - the name of the maintainer (or twitter handle) and a URL (shown in auspice footer)
* `output_folder`: (default: "prepared") will be created inside the current directory and contgain logs + JSONs.
* `segments`: array of strings, or `False` if not segmented...
* `input_format`: (default: fasta) (to do: allow other input formats)
* `input_paths`: {array of strings}. Must be the same number as segments (or 1 if not segmented)
* `header_fields`: dictionary instructing how the fasta header information is encoded. "strain" is essential and "date" is special (see `date_format`). An example:
```
{0:'strain', 2:'accession', 3:'date', 4:'host', 6:'country', 7: 'division'}
```
* `date_format`: (default: `["%Y-%m-%d"]`). Encoding format for date. Can provide multiple formats and they will be tried in order.
* `require_dates` {bool} (default: True) Should sequences without dates be discarded?
* `auspice_filters` {list} (default: []) Those ColorBys which auspice should enable as filters. E.g. "region", "country" etc. Authors are handled separately so don't need to be specified here. These values must also be set as `colors` (see below).

#### filtering settings
* `filters`: {Tuple. Of Tuples. With potentially dictionaries inside. `((a, b), (a, b), ...)`}
  * `a`: name of filter (string)
  * `b`: lambda function _or_ dictionary.

If the filter should be applied to each segment identically, just use a lambda (represented by `b` above).
Some examples:
    * `lambda s: s.attributes['date'] >= datetime(2013,1,1).date()`
    * `s.attributes["host"] not in ["laboratoryderived", "watersample"]`
If the filter is specific to the segment, provide a dictionary of lambdas instead, where the key matches those set in `segments` e.g.
```
{
    "HA": lambda s: len(s.seq)>=1500,
    "NA": lambda s: len(s.seq)>=1200
}
```
Some potential pitfalls:
* By the time filters are applied, the header `strain` has been sanitised, such that the names may have changed.
To avoid problems here, it's easiest to follow the example in `H7N9`:
  * define an array `dropped_strains` before `config`
  * use the filter `("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains])`


#### Complete Genomes & Subsampling
* `ensure_all_segments`: {bool} (Default: `True`). Should only samples with sequences in each segment be included? (If there's only one segment this doesn't have any effect)

#### subsampling
Subsampling is rather complicated in order that most conceivable methods may be employed.
In essence, there are three components:
* categorize the sequences, for instance into year-month groupings (the default)
* assign priorities to each sequence (within each category) - by default this is random
* threshold the sequences (within each category). By default, the 5 sequences with the highest priority in each category are taken.

The `subsample` dict (setting it to `False`, the default, bypasses subsampling) defines three functions to achieve the above three tasks.
These functions can either be lambdas, which are applied to each sequence, or they can be higher order functions which are run when subsampling is started and return a lambda function.

The following keys are searched in the `subsample` dict of the `config` file:
* `category`: can either be a function (usually a lambda) which takes a sequence object and returns a category (anything that can be used as a key in a dict) _OR_ a higher order function that takes `self` (of sequence_set) and returns such a lambda.

* `priority`: similar to above, a lambda with 1 argument (`seq`) _OR_ a higher order function which returns such a lambda. The lambda's return value should be numeric as this is used for sorting.

* `threshold`: an integer (e.g. take _n_ sequences from each category), _OR_ a lambda with 1 argument: the category (see above) of the given sequence, _OR_ a higher order function which returns such a lambda. Lambda's should return an integer.


#### Reference(s)
The Reference sequences is needed for alignment and identification of genes etc, but it doesn't have to be included in the analysis.
The format by which they are defined is rather verbose and hopefully will become part of the database in the future.

**For segmented viruses**:
  * `references`: {dict} with keys corresponding to `segments`. Each value is a {dict} with keys:
    * `path` {string} path to genbank file
    * `include` {int}
    `0`: the reference will be excluded
    `1`: the reference will be added to the pool of sequences, but may be subsampled etc
    `2`: the reference (of that segment) will be included regardless of subsampling.
    Note that currently it may be filtered later on by TreeTime (TODO).
    * `metadata` {dict} of _attribute_ -> _value_, where the attributes are those of `header_fields`.
    `strain` is the only essential attribute, but if you want to add the reference to the dataset then everything should be specified.
    This data (apart from `strain`) is only used if `include > 0` and the strain is not already in the input fasta file.
    * `genes` {`False` | array of {str} | dictionary of {str} key/value pairs} (Default: `False`)
    Genes matching these GenBank annotations will be taken forward for analysis.
    If an array is given, genes will be stored by their GenBank names. If a
    dictionary is given, genes will be mapped from the GenBank names in the keys
    to preferred names in the values.

**For non-segmented viruses:**
  * `reference`: {dict} with keys `path`, `metadata` etc

**Reference Genbank File:**
The reference file must be genbank and only entries with `gene` defined are taken.
This is never as trivial as it should be.
See the references in `H7N9` or `zika` for working examples.


#### colours
Colours are defined for certain fields / traits - normally the same traits that are inferred for nodes in the tree such as `country`, `region`, `host` e.t.c.
Default colour maps will be created for any attributes set here, however you can also supply a file containing custom HEX values.
  * `colors`: False or list of traits (appearing in `header_fields`). If traits are selected, then a section in the JSON will be created with traits -> hex values. By default they will be created using the viridis scale.
  * `color_defs`: _file path_ or [_array_, _of_, _file_, _paths_] - tab separated file(s) joining `trait -> name -> color (hex)`.
  Comment lines start with `#`. See `zika/colors.tsv` for an example.

**To do:** Write a script to change these after auspice JSON creation.

#### latitude & longitude
Similar to colours, these are needed if the data is to be pushed into auspice.
Unlike colours, a file must be provided.
  * `lat_longs`: _False_ or list of traits (appearing in `header_fields`)
  * `lat_long_defs` _file path_ or [_array_, _of_, _file_, _paths_] (default: `'../../fauna/source-data/geo_lat_long.tsv'`)
  ```
  location	country_code	latitude	longitude
  africa	XX	4.070194	21.824559
  north_africa	XX	27.3989987	12.8575109
  subsaharan_africa	XX	0.7603296	25.0743129
  ```
