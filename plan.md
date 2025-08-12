# Test Plan for `subsample.py`

This plan proposes scenarios for two suites:

- **Unit tests (pytest):** focus on Python functions/classes in `subsample.py`.
- **Functional tests (cram):** exercise the `augur subsample` CLI end‑to‑end with small fixtures.

Implementation details and conventions (fixtures layout, helpers, mocks) can follow the wider codebase; below is only the list of scenarios to cover.

---

## Unit tests (pytest)

### 1) `_parse_config()`
- **Parses YAML normally**: returns a dict for a valid config containing `samples`.
- **Treats timestamps as strings**: values that look like dates (e.g., `2020-03-01`) are not auto‑coerced; assert type is `str`.
- **`--config-root` happy path**: when `config_root` is provided and exists, returns that subtree only.
- **`--config-root` missing key**: raises `AugurError` with a helpful message.
- **YAML syntax error**: raises `AugurError` (wrap from `yaml.YAMLError`).
- **Schema validation success**: calls `load_json_schema("schema-subsample-config.json")` and `validate_json(...)` once; no exception.
- **Schema validation failure**: propagates as `AugurError("Config validation failed: …")`.

### 2) `_add_to_args()`
- **Scalar**: appends `[flag, str(value)]`.
- **List/Tuple**: appends `[flag, *map(str, values)]` preserving order.
- **Boolean true with two‑flag tuple**: appends the true flag only.
- **Boolean false with two‑flag tuple**: appends the false flag when present.
- **Boolean false with no false‑flag**: appends nothing.
- **Idempotent ordering**: multiple calls preserve argument order.

### 3) `Sample._construct_filter_args()`
- **Always includes**: `--skip-checks`, `--nthreads 1`, and `--output-strains <tmpfile>`.
- **Extends global args**: all items in `global_filter_args` are present.
- **Maps YAML keys to filter flags**: each key in `SAMPLE_CONFIG` maps to the corresponding flag(s).
- **Boolean mapping**: `probabilistic_sampling=True/False` selects the correct true/false flag.
- **List mapping**: `group_by`, `group_by_weights`, `query_columns` accept list/tuple values.

### 4) `Sample.run()` *
- **Invokes subprocess correctly**: calls `subprocess.run([*augur(shell=False), 'filter', *filter_args], capture_output=True, text=True, check=True)`.
- **Stderr forwarding**: prints each stderr line via `print_err` prefixed with `[<sample name>]` after the process completes.
- **Process error**: on `CalledProcessError`, raises `AugurError` with dedented message containing original stderr and sample name.

### 5) `Sample.cleanup()`
- **Removes temp file**: the temporary `output_strains` file is unlinked.
- **Robust on double‑cleanup**: calling twice doesn’t raise (or is guarded) — decide per convention and assert expected behavior.

### 6) `_run_final_filter()` *
- **Registers and runs**: `augur_filter.register_arguments` receives an `ArgumentParser` and `augur_filter.run` is called with parsed args matching the constructed `filter_args`.

### 7) `run()` orchestration *
- **Builds global args**: only CLI options provided in `args` are included; verify mapping from `GLOBAL_CLI_OPTIONS`.
- **Builds samples**: constructs one `Sample` per config sample and collects their `output_strains` for the final include list.
- **Final args**: begins with `--exclude-all`, then `--include <file...>`, then `--nthreads N`, then global args, then final‑only flags per `FINAL_CLI_OPTIONS`.
- **Parallel execution**: uses `ThreadPoolExecutor(max_workers=args.nthreads)`; all `Sample.run()` futures are awaited.
- **Final filter called**: `_run_final_filter` invoked exactly once after intermediates finish.
- **Cleanup always happens**: `Sample.cleanup()` called for each sample in a `finally`, even on errors.
- **Error bubbling**: if any intermediate sample fails, `run()` raises `AugurError` and still cleans up.

* Note: run functions (`Sample.run()`, `_run_final_filter()`, and `run()`) may be more suitable for functional tests. Feel free to skip these if testing them is unwieldy in pytest.

---

## Functional tests (cram)

These exercise `augur subsample` as a CLI. Use tiny fixture datasets (FASTA + metadata) and minimal configs.

### Happy paths
- **Minimal single‑sample run**: one sample selecting a known set (e.g., `include` IDs) → verify output files are produced; contents match expected union.
- **Multiple samples union**: two samples with overlapping IDs → final output contains the union (no duplicates).
- **Group sampling**: `group_by` + `sequences_per_group` → verify per‑group caps are respected.
- **Deterministic subsampling**: with `--subsample-seed` set, repeated runs yield identical selected IDs.
- **Pass‑through of input options**: using `--sequence-index`, `--metadata-delimiters`, `--metadata-id-columns`, `--metadata-chunk-size` → command succeeds and options reach `augur filter` (observable via debug/log output or side effects per convention).
- **Final outputs**: `--output-sequences` and/or `--output-metadata` are written and non‑empty when appropriate.

### Filtering knobs
- **Include/Exclude by IDs**: `include` and `exclude` lists behave as expected.
- **Where‑clauses**: `include_where` / `exclude_where` select records by metadata conditions.
- **Date bounds**: `min_date` / `max_date` filter by date (ensure timestamp strings are interpreted correctly downstream).
- **Length bounds**: `min_length` / `max_length` filter by sequence length.
- **Ambiguous dates**: `exclude_ambiguous_dates_by` removes records with insufficient date resolution.
- **Non‑nucleotide handling**: `non_nucleotide` toggled (true) to exercise the boolean flag path.
- **Query expression**: `query` with `query_columns` limits evaluation scope.
- **Probabilistic sampling toggle**: run with `probabilistic_sampling` true and false to verify both codepaths.

### Error handling & edge cases
- **Bad YAML**: malformed config → friendly error message from `AugurError`.
- **Unknown `--config-root`**: root key missing → error message identifies the missing key and file.
- **Schema violation**: config missing required keys or with wrong types → schema validation error shown.
- **Intermediate sample failure**: make one sample intentionally fail (e.g., invalid `include_where`) → whole command fails and temp `sample_*` files are cleaned up.
- **Empty selection**: a sample that selects zero IDs still allows overall union to succeed (or yields empty output when all are empty) with a clear message.
- **Missing inputs**: absent `--metadata` and `--sequences` produce a meaningful error from `augur filter` (documented behavior).
