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
These options are described here in sections, but they form one big dictionary - see `H7N9/H7N9.prepare.py` for a working example.
Eventually there will be a "default" config file, so you only have to modify things as necessary.

#### general settings
* `dir`: the current directory - not _augur_ but the virus itself
* `file_prefix`: string used to name the JSONs
* `output_folder`: will be created inside the current directory
* `segments`: array of strings, or `False` if not segmented...
* `input_format`: must be "fasta" currently
* `input_paths`: array of paths (strings). Must be the same number as segments (or 1 if not segmented)
* `header_fields`: dictionary instructing how the fasta header information is encoded. "strain" is essential and "date" is special (see `date_format`). An example:
```
{0:'strain', 2:'accession', 3:'date', 4:'host', 6:'country', 7: 'division'}
```
* `date_format`: encoding format for date - typically `["%Y-%m-%d"]`. Can provide multiple formats and they will be tried in order.
"require_dates": True,

#### filtering settings
* `filters`: Tuple. Of Tuples. With potentially dictionaries inside. `((a, b), (a, b), ...)`
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
* `ensure_all_segments`: bool. Should only samples with sequences in each segement be included?
* `subsample`: `False` or a dict. Dict has keys (each is optional or can have the value `None`)
  * `category`  -- callable that assigns each sequence to a category for subsampling
  * `priority`  -- callable that assigns each sequence a priority to be included in the final sample. this is applied independently in each category
  * `threshold` -- integer (number per category) or callable that determines the number of sequences from each category that is included in the final set. takes arguments, cat and seq.

#### References
References are set here, and stored in the same format as the other sequences, but in `JSON.reference.refName`.
Due to the inconsistencies in genbank / reference formats, you must specify the values for the `header_fields` here.
Perhaps one day it'll be possible to define genes here as well.
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
