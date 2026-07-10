"""
Helpers for YAML-based configuration files.
"""
from argparse import Action, Namespace
from typing import Any, Sequence, TypeAlias, get_args, get_origin

import jsonschema
import yaml

from .argparse_ import ExtendOverwriteDefault
from .errors import AugurError


CONFIG_CLI_OPTIONS_ATTR = "_augur_config_cli_options"
CONFIG_CLI_OPTION_DESTS_ATTR = "_augur_config_cli_option_dests"
ConfigOption: TypeAlias = Any | tuple[Any, str]
ConfigOptionTypes: TypeAlias = dict[str, ConfigOption]


def record_config_cli_option(namespace: Namespace, option_string: str | None, dest: str) -> None:
    """
    Record a config-backed CLI option that was explicitly supplied.

    The original option string is kept for user-facing error messages, while
    the argparse destination is kept for matching the option to YAML keys.
    """
    if option_string is None:
        return

    # Track the exact CLI option string for user-facing conflict errors.
    options = set(getattr(namespace, CONFIG_CLI_OPTIONS_ATTR, set()) or set())
    options.add(option_string)
    setattr(namespace, CONFIG_CLI_OPTIONS_ATTR, options)

    # Map that option string back to its argparse destination so config keys
    # can be compared to CLI options even when their spellings differ.
    option_dests = dict(getattr(namespace, CONFIG_CLI_OPTION_DESTS_ATTR, {}) or {})
    option_dests[option_string] = dest
    setattr(namespace, CONFIG_CLI_OPTION_DESTS_ATTR, option_dests)


class RecordConfigStore(Action):
    """
    Store an argument value and record that the option was explicitly supplied.
    """
    def __call__(
        self,
        parser: Any,
        namespace: Namespace,
        value: Any,
        option_string: str | None = None,
    ) -> None:
        """
        Store the parsed value and remember the option for config conflict checks.
        """
        record_config_cli_option(namespace, option_string, self.dest)
        setattr(namespace, self.dest, value)


class RecordConfigStoreConst(Action):
    """
    Store a constant value and record that the option was explicitly supplied.
    """
    def __init__(
        self,
        option_strings: Sequence[str],
        dest: str,
        const: Any = None,
        default: Any = None,
        required: bool = False,
        help: str | None = None,
        metavar: Any = None,
    ) -> None:
        """
        Initialize a zero-argument action that stores ``const``.
        """
        super().__init__(option_strings, dest, nargs=0, const=const, default=default,
                         required=required, help=help, metavar=metavar)

    def __call__(
        self,
        parser: Any,
        namespace: Namespace,
        value: Any,
        option_string: str | None = None,
    ) -> None:
        """
        Store the constant value and remember the option for config conflict checks.
        """
        record_config_cli_option(namespace, option_string, self.dest)
        setattr(namespace, self.dest, self.const)


class RecordConfigExtendOverwriteDefault(ExtendOverwriteDefault):
    """
    Extend values while recording that the option was explicitly supplied.
    """
    def __call__(
        self,
        parser: Any,
        namespace: Namespace,
        value: Any,
        option_string: str | None = None,
    ) -> None:
        """
        Extend the parsed values and remember the option for config conflict checks.
        """
        record_config_cli_option(namespace, option_string, self.dest)
        super().__call__(parser, namespace, value, option_string)


def apply_config(
    args: Namespace,
    option_types: ConfigOptionTypes,
    command_name: str,
) -> Namespace:
    """
    Apply YAML config values to parsed command arguments.

    Config values are applied only when ``args.config`` is set. CLI options are
    allowed alongside ``--config`` unless the same option also appears in the
    config file.
    """
    if not args.config:
        return args

    cli_options = getattr(args, CONFIG_CLI_OPTIONS_ATTR, set()) or set()
    cli_option_dests = getattr(args, CONFIG_CLI_OPTION_DESTS_ATTR, {}) or {}
    config, config_keys = parse_config(args.config, option_types, command_name)
    conflicting_cli_options = [
        option
        for option in cli_options
        if cli_option_dests.get(option) in config_keys
    ]

    if conflicting_cli_options:
        raise AugurError(
            "--config cannot be used with these CLI options: "
            + ", ".join(sorted(conflicting_cli_options))
        )

    for key, option_type in option_types.items():
        _, dest = config_option_spec(key, option_type)
        value = getattr(config, dest)
        if value is not None:
            setattr(args, dest, value)

    return args


def parse_config(filename: str, option_types: ConfigOptionTypes, command_name: str) -> tuple[Namespace, set[str]]:
    """
    Parse a YAML config file into an argparse namespace and config keys.

    ``option_types`` maps accepted YAML keys to the type used to validate their
    values.
    """
    config_data = load_config(filename, command_name)
    schema = config_json_schema(option_types)

    validate_config(config_data, schema, filename, command_name)

    config = Namespace()
    for option, option_type in option_types.items():
        value_type, dest = config_option_spec(option, option_type)
        value = config_data.get(option) if option in config_data else None

        if value is not None:
            if (value_type is list or get_origin(value_type) is list) and not isinstance(value, list):
                # allow scalars to represent single-item lists
                value = [value]

        setattr(config, dest, value)

    config_keys = {
        config_option_spec(key, option_types[key])[1]
        for key in config_data
    }
    if len(config_keys) != len(config_data):
        raise AugurError(
            f"Invalid {command_name} config {filename!r}: multiple config keys map to the same option"
        )

    return config, config_keys


def load_config(filename: str, command_name: str) -> Any:
    """
    Load a YAML config file.

    YAML timestamps are treated as strings to avoid surprising implicit date
    objects in command configs.
    """
    class CustomLoader(yaml.SafeLoader):
        pass

    def string_constructor(loader: yaml.SafeLoader, node: yaml.Node) -> str:
        return loader.construct_scalar(node)

    CustomLoader.add_constructor("tag:yaml.org,2002:timestamp", string_constructor)

    with open(filename) as handle:
        try:
            config = yaml.load(handle, Loader=CustomLoader)
        except yaml.YAMLError as error:
            raise AugurError(f"Invalid {command_name} config {filename!r}: {error}") from error

    return {} if config is None else config


def config_json_schema(option_types: ConfigOptionTypes) -> Any:
    """
    Return a JSON Schema validator for config options.
    """
    schema = {
        "type": "object",
        "additionalProperties": False,
        "properties": {
            option: json_schema_for_type(option_type)
            for option, option_type in option_types.items()
        },
    }
    Validator = jsonschema.validators.validator_for(schema)
    Validator.check_schema(schema)
    return Validator(schema)


def validate_config(config: Any, schema: Any, filename: str, command_name: str) -> None:
    """
    Validate config data against a JSON Schema validator.
    """
    errors = list(schema.iter_errors(config))
    if not errors:
        return

    error = sorted(errors, key=jsonschema.exceptions.relevance)[0]
    raise AugurError(f"Invalid {command_name} config {filename!r}: {validation_error_message(error)}")


def validation_error_message(error: jsonschema.ValidationError) -> str:
    """
    Return a concise validation error message.
    """
    if error.validator == "additionalProperties":
        additional_property = list(error.instance.keys() - error.schema.get("properties", {}).keys())[0]
        return f"unexpected config option {additional_property!r}"

    return error.message


def config_option_spec(option: str, option_type: ConfigOption) -> tuple[Any, str]:
    """
    Return the value type and destination for a config option.
    """
    if isinstance(option_type, tuple):
        value_type, dest = option_type
        return value_type, dest

    return option_type, option


def json_schema_for_type(value_type: Any) -> dict[str, Any]:
    """
    Return a JSON Schema fragment for a Python type annotation.
    """
    value_type, _ = config_option_spec("", value_type)
    origin = get_origin(value_type)

    import typing
    import types
    if origin is typing.Union or origin is getattr(types, "UnionType", None):
        return {"anyOf": [json_schema_for_type(arg) for arg in get_args(value_type)]}

    if value_type is bool:
        return {"type": "boolean"}
    if value_type is int:
        return {"type": "integer"}
    if value_type is float:
        return {"type": "number"}
    if value_type is str:
        return {"type": "string"}
    if value_type is list or origin is list:
        value_args = get_args(value_type)
        schema: dict[str, Any] = {"type": "array"}
        if value_args:
            schema["items"] = json_schema_for_type(value_args[0])
            return {"anyOf": [schema, schema["items"]]}
        return {"anyOf": [schema, {"type": ["string", "number", "integer", "boolean"]}]}

    raise TypeError(f"unsupported config option type {value_type!r}")
