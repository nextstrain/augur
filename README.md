## Introduction

The nextstrain project is an attempt to make flexible informatic pipelines and visualization tools to track ongoing pathogen evolution as sequence data emerges. The nextstrain project derives from [nextflu](https://github.com/blab/nextflu), which was specific to influenza evolution.

nextstrain is comprised of three components:

* [fauna](https://github.com/nextstrain/fauna): database and IO scripts for sequence and serological data
* [augur](https://github.com/nextstrain/augur): informatic pipelines to conduct inferences from raw data
* [auspice](https://github.com/nextstrain/auspice): web app to visualize resulting inferences

## Augur

*Definition: One held to foretell events by omens.*

Augur is the informatic processing pipeline to track evolution from sequence and serological data.  It is broken into two parts, termed _prepare_ and _process_

#### Prepare:
* loads (and checks) sequence data (from fauna), additional metadata and reference sequences
* Handles segmented viruses
* Filters the data according to user-defined parameters
* subsamples the data if desired
* outputs data (per segment) as JSONs

#### Process:
* Runs independently for each segment / set of inputs
* Takes JSONs (from Prepare)
* aligns sequences
* translates genes
* builds a phylogenetic tree (using RAxML)
* infers timings of internal nodes (using TreeTime)
* infers ancestral states
* infer mutation and clade frequency trajectories through time
* exports intermediate files (alignments, trees) and a JSON bundle for visualization in auspice

## Docs:
* [Prepare](docs/prepare.md)
* [Process](docs/process.md)
* [Format of (auspice) output JSONs](docs/auspice_output.md)


## License and copyright

Copyright 2014-2016 Trevor Bedford and Richard Neher.

Source code to nextstrain is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). nextstrain is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
