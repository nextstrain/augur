## Introduction

The nextstrain project is an attempt to make flexible informatic pipelines and visualization tools to track ongoing pathogen evolution as sequence data emerges. The nextstrain project derives from [nextflu](https://github.com/blab/nextflu), which was specific to influenza evolution.

nextstrain is comprised of three components:

* [fauna](https://github.com/nextstrain/fauna): database and IO scripts for sequence and serological data
* [augur](https://github.com/nextstrain/augur): informatic pipelines to conduct inferences from raw data
* [auspice](https://github.com/nextstrain/auspice): web app to visualize resulting inferences

## Augur

*Definition: One held to foretell events by omens.*

Augur is the informatic processing pipeline to track evolution from sequence and serological data.  It is broken into two parts, termed [_prepare_](docs/prepare.md) and [_process_](docs/process.md), which result in [output JSONs for Auspice](docs/auspice_output.md) (click for documentation).

![flowchart](docs/assets/flow.png)

## How to run:
* see the `README.md` files in the respective pathogen's folder (`flu`, `zika` e.t.c.)
* documentation: [prepare](docs/prepare.md), [process](docs/process.md) and [Auspice JSON format](docs/auspice_output.md)

## License and copyright

Copyright 2014-2016 Trevor Bedford and Richard Neher.

Source code to nextstrain is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). nextstrain is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
