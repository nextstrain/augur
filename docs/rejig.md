## Prepare vs Process
`Prepare` is designed to take inputs and produce a coherent JSON for analysis by `Process`.


## Prepare
Prepare handles pretty much all the input data from fauna / other sources. The idea being that all data wrangling is done here, leaving `Process` free to do proper science. The output is a single JSON (per segment) containing all the information `Process` needs.


| Input        | Status           |
| ------------- | ------------- |
| fauna FASTA      | done |
| references      | to do      |
| CSV metadata | to do      |
| lat/longs | to do      |
| colours | to do      |

### Prepare Config file
See `H7N9/prepare.H7N9.py`. All options are defined in the `config` dict, which is passed to the `Prepare` class. If further customization is needed, just create a new class which inherits from `Prepare` and use that instead.

Most of the options are self explanatory but are listed here for completeness. Eventually there will be a "default" config file, so you only have to modify things as necessary.
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
* `ensure_all_segments`: bool. Should only samples with sequences in each segement be included?
* `filters`: Tuple. Of Tuples. With potentially dictionaries inside. `((a, b), (a, b), ...)`
  * `a`: name of filter (string)
  * `b`: lambda function _or_ dictionary. If the filter should be applied to each segment identically, just use a lambda. some examples:
    * `lambda s: s.attributes['date'] >= datetime(2013,1,1).date()`
    * `s.attributes["host"] not in ["laboratoryderived", "watersample"]`
    * `lambda s: s.id not in [...]`
  * If the function should be specific to the segment, provide a dictionary instead, e.g.
  ```
  {
      "HA": lambda s: len(s.seq)>=1500,
      "NA": lambda s: len(s.seq)>=1200
  }
  ```
* `subsample`: `False` or a dict. Dict has keys (each is optional or can have the value `None`)
  * `category`  -- callable that assigns each sequence to a category for subsampling
  * `priority`  -- callable that assigns each sequence a priority to be included in the final sample. this is applied independently in each category
  * `threshold` -- integer (number per category) or callable that determines the number of sequences from each category that is included in the final set. takes arguments, cat and seq.


## Process
