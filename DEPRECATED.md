# Deprecated

These features are deprecated, which means they are no longer maintained and
will go away in a future major version of Augur. They are currently still
available for backwards compatibility, but should not be used in new code.

## `xopen` major version 1

*Deprecated in version 25.1.0 (July 2024). Removed in version 27.0.0 (December 2024).*

## `augur parse` preference of `name` over `strain` as the sequence ID field

*Deprecated in version 24.2.0 (February 2024). Reordered in version 26.0.0 (September 2024).*

Users who have both 'name' and 'strain' fields in their data, and want to favor using the 'name' field should add the following `augur parse` parameter `--output-id-field 'name'`.

## `augur.utils.read_strains`

*Deprecated in version 24.0.0 (January 2024). Planned for removal March 2024 or after.*

This is part of a [larger effort](https://github.com/nextstrain/augur/issues/1011)
to formalize Augur's Python API.

We recognize the existing usage of this function, so it has been moved to
`augur.io.read_strains`.

## `augur export v1`

*Deprecated in version 22.2.0 (July 2023). Planned for [removal](https://github.com/nextstrain/augur/issues/1266)
January 2024 or after.*

`augur export v2` was introduced in Augur version 6.0.0. Migrate by following
the [official guide](https://docs.nextstrain.org/projects/augur/page/releases/migrating-v5-v6.html).

## `augur ancestral --output`

*Deprecated in version 5.2.0 (December 2019). Removed in version 7.0.0 (April
2020).*

`--output` was replaced by `--output-node-data` to accommodate the addition of
`--output-sequences`.
