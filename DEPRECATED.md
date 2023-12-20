# Deprecated

These features are deprecated, which means they are no longer maintained and
will go away in a future major version of Augur. They are currently still
available for backwards compatibility, but should not be used in new code.

## `augur export v1`

*Deprecated in version 22.2.0 (July 2023). Planned for [removal](https://github.com/nextstrain/augur/issues/1266)
January 2024 or after.*

`augur export v2` was introduced in Augur version 6.0.0. Migrate by following
the [official guide](https://docs.nextstrain.org/projects/augur/page/releases/migrating-v5-v6.html).

## `augur ancestral --output`

*Deprecated in version 5.2.0 (December 2019). Removed in version 7.0.0 (April
2020).*

`--output` was been replaced by `--output-node-data` to accommodate the addition
of `--output-sequences`.
