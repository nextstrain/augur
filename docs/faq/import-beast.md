## Importing BEAST MCC trees into augur

This documentation details how to import BEAST MCC trees using `augur import beast`.
Currently this is most useful for producing auspice-compatable output using `augur export`, however in the future we will provide instructions on how to perform additional analysis using other augur tools.

> BEAST 1 & 2 are extremely versitile tools.
We have tested augur on a number of BEAST runs, using both BEAST & BEAST 2, however there may be issues with your particular run.
Please [get in touch](mailto:hello@nextstrain.org) if you encounter any issues.

### Steps
1. Parse the BEAST tree using `augur import beast`. Most of the options are explained below, but run with `--help` to see them all.
This produces a newick tree file and a node-data JSON file containing BEAST traits as well as (temporal) branch lengths.

2. It may be necessary to modify the format of the traits written to the node-data JSON. These are extracted directly from the BEAST-created annotations in the NEXUS file. For example, if you have encoded location or host as integers, then you should map these back to their true values now.

3. Create an `auspice-config.json` file, which is needed for various display options in auspice.
A template is provided as terminal output from step (1), however there is not enough information in the MCC tree to do this automatically.
Pay particular attention to the color-variable types, which can either be "continuous" or "discrete".

4. Extra metadata can be included here -- either as an additional node-data JSON file, or in TSV format.
Any additional metadata must be both specified in the `auspice-config.json` file and provided to `augur export`.

5. Export auspice-compatible JSONs using `augur export`.
A basic example of what options to supply to this command is provided as terminal output from step (1).


### Taxa naming
The BEAST MCC tree, in NEXUS format, may have `Taxlabels` encoded as integers or strings -- the latter is normally accompanied by a `Translate` block mapping tip names to integers which are used in the actual tree block.


### Calculating the root date
BEAST trees are dated in "time units" from the root date.
Helpfully for us, the sample-date is often encoded in the tip name -- for instance they may follow the format `sample_name|accession|host_name|YYYY-MM-DD` -- and we typically utilise this to calculate the root date.

If the dates are not encoded in the tip names, or no tip-names are used, then you will need to provide the date of the most recent tip (in decimal format) via the `--most-recent-tip-date` argument.

If the dates are provided in the tip names, then we use a regex to extract this (`--tip-date-regex`, the default finds "YYYY-MM-DD" at the end of the tip name). The date format may also be specified if needed via `--tip-date-format`, which is interpreted by the python datetime module,  [see here](https://docs.python.org/3/library/datetime.html#strftime-strptime-behavior) for help with these formats.


### BEAST inferred traits
These are all extracted and stored in the node-data JSON.
The terminal output will list the tratis found, e.g.:
```
Parsed BEAST traits:
name                n(internal) n(terminal)
type                273         274
type_confidence     273         274
height              273         274
height_median       273         272
height_confidence   273         272
posterior           273         0
```

### Examples
```
augur import beast --mcc data/MERS_CoV_mcc.tree --output-tree results/mers.new
    --output-node-data results/beast_data.json
augur export v1 --tree results/mers.new --node-data results/beast_data.json
    --auspice-config config/auspice_config.json
    --output-tree auspice/mers_tree.json --output-meta auspice/mers_meta.json
```

```
augur import beast --mcc data/beast.mcc.nex --output-tree results/mers.new
    --output-node-data results/beast_data.json
    --most-recent-tip-date 2018.43
```

A full [example build can be found here][]. (N.b., this build has been
removed from the current version of the repository; the link here is
to the last version prior to its removal.)

[example build can be found here]: https://github.com/nextstrain/augur/blob/73b8a71103bae3a36b49acdee7c7c75a7e69a751/tests/builds/beast_mers/Snakefile
