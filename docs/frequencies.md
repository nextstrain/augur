# Frequencies

## Aims & Terminology:
The change in the frequency of mutations over time (_mutation frequencies_) and the rise and fall of defined clades (_tree_frequencies_) aid in understanding the evolution of influenza and other pathogens.
These are seperate calculations but share terminology and are conceptually similar.

## Terminology
* **pivots** refer to a series of time slices (coded as a `np.ndarray`), with the frequency calculations being performed at each time point. Note that the keyword argument to `estimate_tree_frequencies` may be a float which is (internally) turned into such an array.

## Mutation frequencies
Mutation frequencies are calculated for each nucleotide mutation and each amino acid mutation in `process.seqs.translations.items()`


#### How to call:
`process.estimate_tree_frequencies()` takes a number of keyword parameters:
* `min_freq` {float} (default: `0.01`)

* `stiffness` {float} (default: `20.0`)

* `inertia` {float} (default: `0.0`)

* `pivots` {int | np.array} (default: `24`)

If a float, then `x` pivots are evenly spaced along the time interval of the samples.
Pivots may instead be spaced by a time slice over the analysis via
`pivots = get_pivots_via_spacing()`, which uses the `pivot_spacing` config parameter.

* `region` {string | tuple | list} (default: `"global"`)

If region is a list, `region[0]` is the name (commonly an acronym) and `region[1]` is a list of regions (or a string for a single region) that is used to reduce the alignment in question by matching to  `seq.attributes['region']`.
If region is a string, the alignment is matched similar to above (both the name and the region refer to this)
If region is `None` or `"global"` then all regions are analysed.

* `include_set` {dict} (default: `{}`)

#### Config parameters
* `estimate_mutation_frequencies` {bool} (default: `False`)
* `pivot_spacing` {float} necessary if you want to calculate pivots through `get_pivots_via_spacing`

#### Side Effects
* `process.mutation_frequencies` {dict}
  * Keys: `(region_name {str}, protein {str})` where proteins are `nuc` and those in `process.seqs.translations`.
  * Values: {dict}
    * Keys: {tuple} `(pos {int}, mutation {str})`
    * Values: {numpy.ndarray, same shape as `pivots`}

* `process.mutation_frequency_confidence` same structure as above.

* `process.mutation_frequency_counts`  {dict}
  * Keys: {str} region name
  * Values: {numpy.ndarray, same shape as `pivots`}
  To check: is this not overwritten by each protein?!?

#### Saving / Restoring
**Saving:**
After each call to `estimate_mutation_frequencies` progress is saved to a pickle in the `processed` directory.


**Restoring:** During the first call to `estimate_mutation_frequencies`, this pickle is loaded (if it exists) and subsequent calculations are skipped where possible (i.e. when the key exists).
The pickle is deemed valid if the names of the sequences are identical.
