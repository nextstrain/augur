## Augur Code Explained

This is a more in-depth explanation of how augur works, detailing the functions to run a pipeline (this guide is currently based off the zika pipeline) and the side-effects of each. It is written for people modifying existing pipelines or creating their own. JSON outputs are detailed [here](../README.md).

### Steps in the pipeline:
[This file](pipeline_steps.md) walks through the functions involved in the [zika pipeline](../zika/zika.py). As such, the main object (class `Process`) is called `zika`. Here are the steps in the pipeline:

* [Initialization](pipeline_steps.md#initialization)
* [Load Data](pipeline_steps.md#load-the-input-data)
* [Filter sequences](pipeline_steps.md#filter-sequences)
* [Subsample Sequences](pipeline_steps.md#subsample-sequences)
* [Align & Translate Sequences](pipeline_steps.md#align-and-translate-sequences)
* [Build the Tree](pipeline_steps.md#build-tree)
* [(optional) Save Progress)](pipeline_steps.md#save-progress)
* [Molecular Clock Filtering](pipeline_steps.md#molecular-clock-filtering)
* [Annotate Tree](pipeline_steps.md#annotate-tree)
* [Geographical Inference](pipeline_steps.md#geographical-inference)
* [Export Data](pipeline_steps.md#export-data)


### Input files:
Augur is designed to take files from [fauna](https://github.com/nextstrain/fauna), namely:
* A multi fasta file, with information encoded in the header.
* A genbank reference sequence, which doesn't have to appear in the fasta file.
* (optional) A tsv file used to map regions to lat/long values. Header: `location	country_code	latitude	longitude`
* (optional) a file linking sequence names to additional attributes (used to decorate the tree)

### Organization of augur code
* The main processing file (e.g. `zika/zika.py`) contains the project-specific settings and a number of commands to run the pipeline. A number of augur and treetime classes are used (depending on what's run) most notably:

* The _process_ class (`base/process.py`) contains higher-level functions, which mostly wrap methods from the following classes. Most of the commands in the main processing file are methods from this class.

* The *sequence_set* class (`base/sequences.py`) deals with sequence-oriented functions, such as parsing, aligning, subsampling, translating, filtering etc

* The _tree_ class (`base/tree.py`) builds the tree, performs geo-inference, and interacts with the _TreeTime_ class from [treetime](https://github.com/neherlab/treetime).

* The _titers_ class (`base/titer_model.py`)

* The *tree_predictor* class (`base/prediction.py`)


### Status of these docs:
* Based off augur commit `14549d4`
* To do:
  * Titers
  * frequencies
