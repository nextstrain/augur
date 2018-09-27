# __NEXT__


# 3.0.2.dev1 (27 September 2018)

## Bug fixes

* translate: Fix broken `--help` message


# 3.0.1.dev1 (27 September 2018)

## Features

* align and tree: The --nthreads option now accepts the special value "auto" to
  automatically set the number of threads to the number of CPU cores available.

* Alias `augur --version` to `augur version`

## Bug fixes

* tree: The --nthreads option is now respected.  Previously all tree builders
  were ignoring the value and using either 2 threads (RAxML, IQ-TREE) or as
  many threads as cores (FastTree, if the OpenMP version).

* translate: Check for and, if necessary pad, nucleotide sequences which aren't
  a multiple of 3 earlier to avoid errors later.

* export: Optionally write inferred nucleotide and amino acid sequences (or
  mutations) to a separate file.

* export: Omit genes with no amino acid mutations.

* validate: Allow underscores in gene names.

* refine: Remove unused --nthreads argument.

* ancestral, filter, tree, refine: Exit 1 instead of -1 on error.

* Print the help message, instead of throwing an exception, when `augur` is run
  without arguments.

## Documentation

* Briefly describe each command in its `--help` output and in the global `augur
  --help` output.

* Revamp README to emphasize new, modular augur and make it suitable for
  inclusion on PyPi.

* Reconciled conflicting license declarations; augur is AGPLv3 (not MIT)
  licensed like the rest of Nextstrain.

* Include URLs for bug reports, the change log, and the source on PyPi.

## Data

* Geographic coordinates added for the Netherlands and the Philippines.

## Development

* Reset the `release` branch when rewinding a failed local release process.

* Refactor the augur program and command architecture for improved
  maintainability.


# 3.0.0.dev3 (4 September 2018)

## Development

* Use an allowed Topic classifier so we can upload to PyPi

* Ignore distribution egg-info build files


# 3.0.0.dev2 (4 September 2018)

## Features

* Export: Add safety checks for optional annotations and geo data

* Include more lat/longs in the default geo data

## Development

* Add release tooling

* Document the release process and a few development practices

* Travis CI: Switch to rebuilding the Docker image only for new releases

* Remove ebola, lassa, tb, WNV, and zika builds now in their own repos.  These
  builds are now available at URLs like <https://github.com/nextstrain/ebola>,
  for example.


# 3.0.0.dev1 (unreleased)

## Development

* Start versioning augur beginning with 3.0.0.  A new `augur version` command
  reports the running version.
