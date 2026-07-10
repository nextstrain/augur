# Command-Line Configuration File Specification

## 1. Objective
Enable command-line tools to accept input arguments via a YAML configuration file in addition to standard command-line interface (CLI) flags. This improves flexibility when used in Snakemake workflows.

## 2. Core Design Principles
The design relies on a strict separation of inputs and an explicit mapping between configuration files and the CLI. To achieve this, any implementation must follow these structural rules:

### 2.1 Configuration Mapping
Every command that supports a configuration file must define a static "Configuration Mapping". This mapping serves as the single source of truth connecting YAML options to the CLI.

The mapping must dictate:
1. **Key Name:** The exact string key expected in the YAML file (e.g., `max_iter`).
2. **Type:** The expected data type of the value (e.g., `integer`, `boolean`, `list[str]`).
3. **CLI Equivalent:** The corresponding CLI flag(s) that this YAML key replaces (e.g., `--max-iter`).

**Advanced Mapping Rules:**
- The mapping should simplify logical inversions into single, 3-state boolean variables (e.g., `covariance: True / False` for `--covariance`, `--no-covariance`). A missing key (None) implies the default behavior.

### 2.2 CLI vs. Config State Tracking
To merge configurations safely, the underlying CLI parser must be able to distinguish between an argument that was *explicitly provided* by the user on the CLI versus an argument that simply holds a *default fallback value*. 

The system must track the source of every argument (CLI, Config, or Default).

## 3. Merging and Validation Rules
When a configuration file is passed to a command, the system reads the file and merges it with the CLI arguments according to the following strict rules:

### 3.1 Mutual Exclusion (No Overrides)
To keep invocations simple and predictable, a user may provide an option via the CLI *or* via the configuration file, **but not both**. 
- If an option is specified in the configuration file, and the equivalent flag is also passed explicitly on the CLI, the system must **raise a fatal error**.
- The system does not support "CLI overriding config" or "config overriding CLI".

### 3.2 Strict Validation
- **Unknown Keys:** If the YAML configuration file contains any keys that do not exist in the defined Configuration Mapping, the system must **raise a fatal error**. This protects against typos and silent failures.
- **Type Checking:** The values provided in the YAML file should be validated against the expected types defined in the mapping.

## 4. Automatic Documentation
Because the Configuration Mapping acts as a single source of truth, documentation describing how YAML config keys map to CLI options should be generated automatically. Implementations should provide a mechanism (such as a documentation generation macro or directive) that consumes the mapping and outputs a reference table. This guarantees that documentation never falls out of sync with the supported configuration options.
