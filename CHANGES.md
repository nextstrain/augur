# __NEXT__


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
