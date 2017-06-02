# Prepare
Prepare handles pretty much all the input data from fauna / other sources.
The idea being that all data wrangling is done here, leaving `Process` free to do proper science.
The output is a single JSON (per segment) containing all the information `Process` needs.
The hope is that just modifying the `config` dictionary (in, e.g., `H7N9.prepare.py`) is all that is needed for any analysis.
If you need something a bit more bespoke that can't be done through the config, the idea would be to create a new class which inherits from `Prepare` and use that instead.

### Status:

| Input        | Status           |
| ------------- | ------------- |
| fauna FASTA      | DONE |
| references    | DONE      |
| CSV metadata | to do      |
| lat/longs | DONE      |
| colours (auto) | DONE      |
| colours (user defined) | DONE      |

### `Config` dictionary
As mentioned above, all options are defined in the `config` dict, which is passed to the `Prepare` class.
Many options have defaults provided, and if so are not required to be in the config file.

#### general settings
* `dir`: the current directory - not _augur_ but the virus itself
* `file_prefix`: string used to name the JSONs
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

* `threshold`: an integer (e.g. take _n_ sequences from each category), _OR_ a lambda with 1 argument: they category (see above) of the given sequence, _OR_ a higher order function which returns such a lambda. Lambda's should return an integer.


#### References
References are set here, and stored in the same format as the other sequences, but in `JSON.reference.refName`.
Due to the inconsistencies in genbank / reference formats, you must specify the values for the `header_fields` here.
Perhaps one day it'll be possible to define genes here as well.
Ideally these will be provided in the database so this will be unnecessary.
It's not essential to have a reference for each segment, but it's highly recommended.
  * `references`: Dictionary with keys corresponding to `segments`. Each value is a dictionary with:
    * `path`: path to genbank file
    * `metadata`: dictionary of `country`->`XXX` etc. Corresponds to the values in `header_fields`. Missing data is given "unknown"
    * `use`: bool. The reference is used in the alignment stage to define coding regions. But it's not required to be used for the phylogeny etc, and in some cases it's useful not to use it here (e.g. it's older than the samples in question, too divergent e.t.c). Default: False
    * `genes`: False or array of strings. Genes matching these will be taken forward for analysis. Default: False

If it's a non-segmented virus, use a single dictionary at `reference` (not `references`)

The reference file must be genbank and only entries with `gene` defined are taken. This is never as trivial as it should be. See the references in `H7N9` for working examples.


#### colours
Colours (sic) are defined for certain fields / traits - usually `country` and `region`
These are crucial if you want to push your data into auspice.
Eventually i'll write a script to change these after the JSONs have been created.
But for now...
  * `colors`: False or list of traits (appearing in `header_fields`). If traits are selected, then a section in the JSON will be created with traits -> hex values. By default they will be created using the viridis scale.
  * `color_defs`: _file path_ or [_array_, _of_, _file_, _paths_] - tab separated file(s) joining trait -> color (hex), or dictionary linking trait values to hexes

#### latitude & longitude
Similar to colours, these are needed if the data is to be pushed into auspice.
Unlike colours, a file must be provided.
  * `lat_longs`: _False_ or list of traits (appearing in `header_fields`)
  * `lat_long_defs`: _file path_ or [_array_, _of_, _file_, _paths_]
  ```
  location	country_code	latitude	longitude
  africa	XX	4.070194	21.824559
  north_africa	XX	27.3989987	12.8575109
  subsaharan_africa	XX	0.7603296	25.0743129
  ```
