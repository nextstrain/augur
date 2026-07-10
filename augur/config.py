"""
Support for YAML configuration files that provide CLI argument values.

Commands opt in individual arguments by using :func:`add_config_argument`
instead of ``parser.add_argument``.  At the end of ``register_parser``, call
:func:`apply_sentinels` to finalize the mapping.  In ``run()``, call
:func:`merge_config` to load, validate, and merge the config file.
"""
import yaml
from dataclasses import dataclass, field
from typing import Any, List

from .errors import AugurError


class _UnsetType:
    """Sentinel indicating an argument was not explicitly provided on the CLI."""
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __repr__(self):
        return '_UNSET'

    def __bool__(self):
        return False

_UNSET = _UnsetType()


@dataclass
class ConfigOption:
    """A CLI argument that can also be set via a YAML config file."""
    key: str
    dest: str
    cli_flags: List[str]
    type: str  # "int", "float", "bool", "str", "list[int]", "list[str]"
    default: Any = field(default=_UNSET, repr=False)


def add_config_argument(group, *args, config_key=None, **kwargs):
    """Wrapper around ``group.add_argument`` that marks the argument as
    config-backed.

    *config_key* is the expected YAML key name.  When omitted it defaults to
    ``action.dest``.
    """
    action = group.add_argument(*args, **kwargs)
    action.config_key = config_key if config_key is not None else action.dest
    return action


def apply_sentinels(parser):
    """Scan *parser* for config-backed arguments, build the mapping, swap
    defaults with :data:`_UNSET`, and store the mapping on *parser* via
    ``set_defaults``.

    Must be called **after** all ``add_argument`` and ``set_defaults`` calls.
    """
    mapping = []
    seen_keys = {}  # config_key -> ConfigOption

    for action in parser._actions:
        config_key = getattr(action, 'config_key', None)
        if config_key is None:
            continue

        if config_key in seen_keys:
            # Merge flags for shared-dest entries (e.g. --covariance / --no-covariance)
            seen_keys[config_key].cli_flags.extend(action.option_strings)
            continue

        entry = ConfigOption(
            key=config_key,
            dest=action.dest,
            cli_flags=list(action.option_strings),
            type=_infer_type(action),
        )
        mapping.append(entry)
        seen_keys[config_key] = entry

    # Capture original defaults, then replace with sentinel
    for entry in mapping:
        entry.default = parser.get_default(entry.dest)
        parser.set_defaults(**{entry.dest: _UNSET})

    parser.set_defaults(_config_mapping=mapping)


def merge_config(args, command_name):
    """Load a YAML config (if provided), detect conflicts, merge values, and
    restore defaults for any arguments that were set by neither the CLI nor the
    config.

    Call this at the top of ``run(args)``::

        merge_config(args, 'refine')
    """
    mapping = args._config_mapping
    config_path = getattr(args, 'config', None)

    config = {}
    if config_path:
        config = _load_config(config_path, mapping, command_name)

        # Conflict detection: config key present AND CLI arg explicitly set
        conflicts = []
        for entry in mapping:
            if entry.key in config and getattr(args, entry.dest) is not _UNSET:
                conflicts.append(entry.cli_flags[0])
        if conflicts:
            raise AugurError(
                f"--config cannot be used with these CLI options: {', '.join(conflicts)}"
            )

        # Apply config values
        for entry in mapping:
            if entry.key in config:
                setattr(args, entry.dest, config[entry.key])

    # Restore original defaults for anything still unset
    for entry in mapping:
        if getattr(args, entry.dest) is _UNSET:
            setattr(args, entry.dest, entry.default)


def _load_config(path, mapping, command_name):
    """Read a YAML config file and validate its keys and value types."""
    try:
        with open(path) as f:
            raw = yaml.safe_load(f)
    except FileNotFoundError:
        raise AugurError(f"Invalid {command_name} config {path!r}: file not found")
    except yaml.YAMLError as e:
        raise AugurError(f"Invalid {command_name} config {path!r}: {e}")

    if raw is None:
        return {}

    if not isinstance(raw, dict):
        raise AugurError(
            f"Invalid {command_name} config {path!r}: expected a YAML mapping, got {type(raw).__name__}"
        )

    valid_keys = {entry.key for entry in mapping}
    for key in raw:
        if key not in valid_keys:
            raise AugurError(
                f"Invalid {command_name} config {path!r}: unexpected config option {key!r}"
            )

    config = {}
    entry_by_key = {entry.key: entry for entry in mapping}
    for key, value in raw.items():
        entry = entry_by_key[key]
        config[key] = _validate_value(key, value, entry.type, path, command_name)

    return config


def _validate_value(key, value, expected_type, path, command_name):
    """Validate and coerce a single config value against *expected_type*."""
    prefix = f"Invalid {command_name} config {path!r}: option {key!r}"

    if expected_type == 'bool':
        if not isinstance(value, bool):
            raise AugurError(f"{prefix} expects a boolean, got {type(value).__name__}")
        return value

    if expected_type == 'int':
        if not isinstance(value, int) or isinstance(value, bool):
            raise AugurError(f"{prefix} expects an integer, got {type(value).__name__}")
        return value

    if expected_type == 'float':
        if isinstance(value, bool):
            raise AugurError(f"{prefix} expects a number, got bool")
        if not isinstance(value, (int, float)):
            raise AugurError(f"{prefix} expects a number, got {type(value).__name__}")
        return float(value)

    if expected_type == 'str':
        # Accept scalars and convert to string (e.g. coalescent: 0.01)
        if isinstance(value, bool):
            raise AugurError(f"{prefix} expects a string, got bool")
        if not isinstance(value, (str, int, float)):
            raise AugurError(f"{prefix} expects a string, got {type(value).__name__}")
        return str(value)

    if expected_type in ('list[str]', 'list[int]'):
        element_type = str if expected_type == 'list[str]' else int

        # Auto-wrap bare scalar into single-element list
        if not isinstance(value, list):
            value = [value]

        for i, item in enumerate(value):
            if element_type is str:
                if isinstance(item, bool):
                    raise AugurError(f"{prefix} expects a list of strings, got bool at index {i}")
                if not isinstance(item, (str, int, float)):
                    raise AugurError(f"{prefix} expects a list of strings, got {type(item).__name__} at index {i}")
                value[i] = str(item)
            else:
                if not isinstance(item, int) or isinstance(item, bool):
                    raise AugurError(f"{prefix} expects a list of integers, got {type(item).__name__} at index {i}")
        return value

    raise ValueError(f"Unknown config type {expected_type!r}")


def _infer_type(action):
    """Infer the config value type from an argparse action's properties."""
    # Boolean store actions
    if action.const is True or action.const is False:
        return 'bool'

    # Scalar type from action.type
    if action.type is int:
        scalar = 'int'
    elif action.type is float:
        scalar = 'float'
    else:
        scalar = 'str'

    # List types
    if action.nargs in ('+', '*'):
        return f'list[{scalar}]'

    return scalar
