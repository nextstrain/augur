# Development Docs for Contributors

Thank you for helping us to improve Nextstrain! This document describes:

- Getting Started
- Contributing code
  - Running local code changes
  - Testing
  - Creating a release
  - Continuous integration
- Contributing documentation
  - Formats (Markdown and reStructuredText)
  - Documentation structure
  - Building documentation

## Getting started

To be an effective, productive contributor, please start by reading the
[**Nextstrain contributing guide**](https://github.com/nextstrain/.github/blob/master/CONTRIBUTING.md)
for useful information about how to pick an issue, submit your contributions, and so on.

This project strictly adheres to the
[Contributor Covenant Code of Conduct](https://github.com/nextstrain/.github/blob/master/CODE_OF_CONDUCT.md).

Please see the [project board](https://github.com/orgs/nextstrain/projects/6) for currently available issues.

## Contributing code

We currently target compatibility with Python 3.6 and higher. As Python releases new versions,
the minimum target compatibility may be increased in the future.

Versions for this project, Augur, from 3.0.0 onwards aim to follow the
[Semantic Versioning rules](https://semver.org).

### Running local changes

While you are making code changes, you will want to run augur to see it behavior with those changes.
To to test your local changes (without installing them to your system), run the following
**convenience script** from the root of your cloned git repository:

```bash
./bin/augur
```

Note that the `./bin/augur` convenience script is not installing `augur` system-wide with pip.

As an alternative to using the convenience script, you can install augur from source
as an **editable package** so that your global `augur` command always uses your
local source code copy:

```bash
pip install -e '.[dev]'
```

Using an "editable package" is not recommended if you want to be able to compare output
from a stable, released version of augur with your development version (e.g. comparing
output of `augur` installed with pip and `./bin/augur` from your local source code).

### Testing

Writing good tests and running tests helps maintain code quality and eases future refactoring.
We use [pytest](https://docs.pytest.org) and [Cram](https://bitheap.org/cram/) to test augur.
This section will describe briefly:

- Writing tests
  - Unit tests
  - Doctests
  - Functional tests
- Running tests
  - Locally
  - Continuous Integration

#### Writing Tests

It's good practice to write **unit tests** for any code contribution.
The [pytest documentation](https://docs.pytest.org) and [Python documentation](https://docs.python.org) are good references for unit tests.
Augur's unit tests are located in the `tests` directory and there is generally one test file for each code file.

On the other hand, [**doctests**](https://docs.python.org/3/library/doctest.html) are a type of tests that are written within a module's docstrings.
They can be helpful for testing a real-world example and determining if a regression is introduced in a particular module.

A pull request that contributes new code **should always contain unit tests**.
Optionally, a pull request may also contain doctests if the contributor believes a doctest would improve the documentation and execution of a real world example.

We test augur's command line interface with functional tests implemented with the [Cram framework](https://bitheap.org/cram/).
These tests complement existing unit tests of individual augur Python functions by running augur commands on the shell and confirming that these commands:

1. execute without any errors
2. produce exactly the expected outputs for the given inputs

These tests can reveal bugs resulting from untested internal functions or untested combinations fo internal functions.

Functional tests should either:

* suitably test a single augur command with an eponymously named Cram file in `tests/functional/` (e.g., `mask.t` for augur mask)

OR

* test a complete build with augur commands with an appropriately named Cram file in `tests/builds/` (e.g., `zika.t` for the example Zika build)

##### Functional tests of specific commands

Functional tests of specific commands consist of a single Cram file per test and a corresponding directory of expected inputs and outputs to use for comparison of test results.

The Cram file should test most reasonable combinations of command arguments and flags.

##### Functional tests of example builds

Functional tests of example builds use output from a real Snakemake workflow as expected inputs and outputs.
These tests should confirm that all steps of a workflow can execute and produce the expected output.
These tests reflect actual augur usage in workflows and are not intended to comprehensively test interfaces for specific augur commands.

The Cram file should replicate the example workflow from start to end.
These tests should use the output of the Snakemake workflow (e.g., files in `zika/results/` for the Zika build test) as the expected inputs and outputs.

##### Comparing outputs of augur commands

Compare deterministic outputs of augur commands with a `diff` between the expected and observed output files.
For extremely simple deterministic outputs, use the expected text written to standard output instead of creating a separate expected output file.

To compare trees with stochastic branch lengths:

1. provide a fixed random seed to the tree builder executable (e.g., `--tree-builder-args "-seed 314159"` for the “iqtree” method of augur tree)
2. use `scripts/diff_trees.py` instead of `diff` and optionally provide a specific number to  `--significant-digits` to limit the precision that should be considered in the diff

To compare JSON outputs with stochastic numerical values, use `scripts/diff_jsons.py` with the appropriate `--significant-digits` argument.

Both tree and JSON comparison scripts rely on [deepdiff](https://deepdiff.readthedocs.io/en/latest/) for underlying comparisons.

#### Running Tests

You've written tests and now you want to run them to see if they are passing.
First, you will need to [install the complete Nextstrain environment](https://nextstrain.org/docs/getting-started/local-installation) and augur dev dependencies as described above.
Next, run all augur tests with the following command from the root, top-level of the augur repository:

```bash
./run_tests.sh
```

For rapid execution of a subset of unit tests (as during test-driven development), the `-k` argument will disable code coverage and functional tests and pass directly to pytest to limit the tests that are run.
For example, the following command only runs unit tests related to augur mask.

```bash
./run_tests.sh -k test_mask
```

Troubleshooting tip: As tests run on the development code in the augur repository, your environment should not have an existing augur installation that could cause a conflict in pytest.

We use continuous integration with Travis CI to run tests on every pull request submitted to the project.
We use [codecov](https://codecov.io/) to automatically produce test coverage for new contributions and the project as a whole.

### Releasing

New releases are tagged in git using an "annotated" tag.  If the git option
`user.signingKey` is set, the tag will also be [signed][].  Signed tags are
preferred, but it can be hard to setup GPG correctly.  The `release` branch
should always point to the latest release tag.  Source and wheel (binary)
distributions are uploaded to [the nextstrain-augur project on
PyPi](https://pypi.org/project/nextstrain-augur).

There is a `./devel/release` script which will prepare a new release from your
local repository.  It ends with instructions for you on how to push the release
commit/tag/branch and how to upload the built distributions to PyPi.  You'll
need [a PyPi account][] and [twine][] installed to do the latter.

[signed]: https://git-scm.com/book/en/v2/Git-Tools-Signing-Your-Work
[a PyPi account]: https://pypi.org/account/register/
[twine]: https://pypi.org/project/twine

### Travis CI

Branches and PRs are tested by Travis CI jobs configured in `.travis.yml`.

Our Travis config uses two build stages: _test_ and _deploy_.  Jobs in the
_test_ stage always run, but _deploy_ jobs only run sometimes (see below).

The set of _test_ jobs are explicitly defined instead of auto-expanded from the
implicit job property matrix. Since top-level properties are inherited by all
jobs regardless of build stage, making the matrix explicit is less confusing
and easier to reason about. YAML's anchor (`&foo`) and alias merge key (`<<:
*foo`) syntax let us do this without repeating ourselves unnecessarily.

New releases, via pushes to the `release` branch, trigger a new [docker-base][]
build to keep the Docker image up-to-date. This trigger is implemented in the
_deploy_ stage, which is implicitly conditioned on the previous _test_ stage's
successful completion and explicitly conditioned on a non-PR trigger on the
`release` branch. Note that currently we cannot test this _deploy_ stage
without making a release.

It can sometimes be useful to verify the config is parsed as you expect using
<https://config.travis-ci.com/explore>.

[docker-base]: https://github.com/nextstrain/docker-base

## Contributing documentation

[Documentation](https://nextstrain-augur.readthedocs.io) is built using [Sphinx](http://sphinx-doc.org/) and hosted on [Read The Docs](https://readthedocs.org/).
Versions of the documentation for each augur release and git branch are available and preserved.
Read The Docs is updated automatically from commits and releases on GitHub.

### Doc Formats

Documentation is mostly written as [reStructuredText](http://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html) (.rst) files, but they can also be Markdown (.md) files.
There are advantages to both formats:

- reStructuredText enables python-generated text to fill your documentation as in the
  auto-importing of modules or usage of plugins like `sphinx-argparse` (see below).
- Markdown is more intuitive to write and is widely used outside of python development.
- If you don't need autogeneration of help documentation, then you may want to stick with writing Markdown.

Sphinx, coupled with reStructuredText, can be tricky to learn.
Here's a [subset of reStructuredText worth committing to memory](https://simonwillison.net/2018/Aug/25/restructuredtext/) to help you get started writing these files.

Many Sphinx reStructuredText files contain a [directive](http://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html) to add relations between single files in the documentation known as a Table of Contents Tree ([TOC Tree](https://documentation.help/Sphinx/toctree.html)).

Human-readable augur and augur subcommand documentation is written using a Sphinx extension called [sphinx-argparse](https://sphinx-argparse.readthedocs.io/en/latest/index.html).

### Folder structure

The documentation source-files are located in `./docs`, with `./docs/index.rst` being the main entry point.
Each subsection of the documentation is a subdirectory inside `./docs`.
For instance, the tutorials are all found in `./docs/tutorials` and are included in the documentation website via the directive in `./docs/index.rst`.

### Building documentation

Building the documentation locally is useful to test changes.
First, make sure you have the development dependencies of augur installed:

```bash
pip install -e '.[dev]'
```

This installs packages listed in the `dev` section of `extras_require` in _setup.py_,
as well as augur's dependencies as necessary.

Sphinx and make are used when **building** documentation. Here are some examples that you
may find useful:

Build the HTML output format by running:

```bash
make -C docs html
```

Sphinx can build other formats, such as epub. To see other available formats, run:

```bash
make -C docs help
```

To update the API documentation after adding or removing an augur submodule,
autogenerate a new API file as follows.

```bash
sphinx-apidoc -T -f -MeT -o docs/api augur
```

To make doc rebuilds faster, Sphinx caches built documentation by default,
which is generally great, but can cause the sidebar of pages to be stale.
You can clean out the cache with:

```bash
make -C docs clean
```

To **view** the generated documentation in your browser, Mac users should run:

```bash
open docs/_build/html/index.html
```

Linux users can **view** the docs by running:

```bash
xdg-open docs/_build/html/index.html
```

This will open your browser for you to see and read your work.
