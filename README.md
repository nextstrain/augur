## Introduction

Nextstrain is an open-source project to harness the scientific and public health potential of pathogen genome data. We provide a continually-updated view of publicly available data with powerful analytics and visualizations showing pathogen evolution and epidemic spread. Our goal is to aid epidemiological understanding and improve outbreak response.

Nextstrain is comprised of three primary components:

* [fauna](https://github.com/nextstrain/fauna): database and IO scripts for sequence and serological data
* [augur](https://github.com/nextstrain/augur): informatic pipelines to conduct inferences from raw data
* [auspice](https://github.com/nextstrain/auspice): web app to visualize resulting inferences

Resulting data and inferences are available live at the website [nextstrain.org](http://nextstrain.org).

## Augur

*Definition: One held to foretell events by omens.*

Augur is the informatic processing pipeline to track evolution from sequence and serological data.  It is broken into two parts, termed [_prepare_](docs/prepare.md) and [_process_](docs/process.md), which result in [output JSONs for Auspice](docs/auspice_output.md).

![flowchart](docs/assets/flow.png)

## Documentation

[Docs](docs/) are available for [prepare](docs/prepare.md), [process](docs/process.md) and [Auspice JSON format](docs/auspice_output.md).

## Virus builds

Each virus build consists of a `prepare.py` and `process.py` file. Current builds provided are:

* [dengue virus](dengue/)
* [Ebola virus](ebola/)
* [seasonal influenza virus](flu/)
* [avian influenza H7N9](avian/)
* [mumps virus](mumps/)
* [Zika virus](zika/)

## License and copyright

Copyright 2014-2017 Trevor Bedford and Richard Neher.

Source code to Nextstrain is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). Nextstrain is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
