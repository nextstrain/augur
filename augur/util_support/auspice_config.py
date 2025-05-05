import sys
import os
import json
from ..io.file import open_file
from ..errors import AugurError
from .warnings import deprecated
from ..types import ValidationMode
from ..validate import auspice_config_v2 as validate_auspice_config_v2, ValidateError, validation_failure
from collections import defaultdict
from typing import Union, Any
from collections.abc import Callable


# Deprecated attr keys that are automatically updated during `augur export`
# old_key: new_key
DEPRECATED_KEYS = {
    "authors": "author",
    "numdate": "num_date"
}

def update_deprecated_names(name):
    # correct deprecated keys
    return DEPRECATED_KEYS.get(name, name)

def read_json(fname):
    if not (fname and os.path.isfile(fname)):
        print("ERROR: config file %s not found."%fname)
        return defaultdict(dict)

    try:
        with open_file(fname, 'rb') as ifile:
            config = json.load(ifile)
    except json.decoder.JSONDecodeError as err:
        print("FATAL ERROR:")
        print("\tCouldn't parse the JSON file {}".format(fname))
        print("\tError message: '{}'".format(err.msg))
        print("\tLine number: '{}'".format(err.lineno))
        print("\tColumn number: '{}'".format(err.colno))
        print("\tYou must correct this file in order to proceed.")
        sys.exit(2)

    return config

def _merge_lists(base: dict[str, Any], overlay: dict[str, Any], config_key: str, identity_func: Callable) -> dict[str, Any]:
    """
    This merging approach is used for config keys which are encoded as lists
    such as maintainers, colorings, filters.

    Returns an extended list with elements from the base (config's) list and the
    overlay (config's) list. If elements in the overlay are present in the base
    list, as determined by the *identity_func* then the base element is replaced
    with the overlay element.

    NOTE: Currently there is no way to remove elements from the base list, nor
    is it possible to insert a new (overlay) element at a specific position.
    Both are possible but increase the conceptual and programmatic complexity.
    """

    if config_key not in overlay:
        return base
    if config_key not in base:
        return {**base, config_key: overlay[config_key]}

    if not isinstance(overlay[config_key], list) or not isinstance(base[config_key], list):
        raise AugurError(f"Config merging for {config_key!r} failed as one (or more) entries were not lists")

    values = [*base[config_key]]
    for overlay_el in overlay[config_key]:
        for (idx, base_el) in enumerate(values):
            if identity_func(overlay_el, base_el):
                values[idx] = overlay_el # replace existing element with one from the overlay config
                break
        else:
            values.append(overlay_el)

    return {**base, config_key: values}

def _merge_scalar(base: dict[str, Any], overlay: dict[str, Any], config_key: str) -> dict[str, Any]:
    return {**base, config_key: overlay[config_key]} if config_key in overlay else base

def _merge_dicts(base: dict[str, Any], overlay: dict[str, Any], config_key: str) -> dict[str, Any]:
    """
    Dicts are _not_ recursively merged, they are simply merged by adding new keys from the overlay to
    those keys already present in the base. Identical keys replace the original key.
    """
    if config_key not in overlay:
        return base
    if config_key not in base:
        return {**base, config_key: overlay[config_key]}

    if not isinstance(overlay[config_key], dict) or not isinstance(base[config_key], dict):
        raise AugurError(f"Config merging for {config_key!r} failed as one (or more) entries were not dictionaries")

    return {**base, config_key: {**base[config_key], **overlay[config_key]}}

def _geo_resolution_id(el: Union[str, dict[str,str]]) -> str:
    return el['key'] if isinstance(el, dict) else el

def _replace_deprecated(config: dict[str,Any], old_name: str, new_name: Union[None,str], modify: Union[None, Callable]=None):
    if old_name not in config:
        return # NO-OP
    elif new_name is None:
        deprecated(f"[config file] key {old_name!r} is no longer used and has been dropped from your config")
        del config[old_name]
    elif old_name == new_name: # indicates deprecations within an unchanged top-level key
        assert modify is not None, "[internal error] modify must be a function"
        config[new_name] = modify(config[old_name])
    elif new_name in config:
        deprecated(f"[config file] top level keys {new_name!r} and {old_name!r} were both present, " \
                   "however the latter is a deprecated version of the former and will be ignored.")
        del config[old_name]
    else:
        deprecated(f"[config file] top level key {old_name!r} is deprecated. Renaming to {new_name!r}" +
                   (" (with modifications)" if modify else ""))
        config[new_name] = modify(config[old_name]) if modify else config[old_name]
        del config[old_name]

def _parse_color_options(options: Any) -> list[dict]:
    # parse v1-style 'color_options' & convert to v2-style 'colorings'
    if not isinstance(options, dict):
        raise AugurError(f"[config file] The deprecated field 'color_options' must be a dictionary")
    colorings = []
    for key, info in options.items():
        # note: if both legendTitle & menuItem are present then we use the latter. See https://github.com/nextstrain/auspice/issues/730
        if "menuItem" in info:
            deprecated("[config file] coloring '{}': 'menuItem' has been replaced with 'title'. Using 'menuItem' as 'title'.".format(key))
            info["title"] = info["menuItem"]
            del info["menuItem"]
        if "legendTitle" in info:
            if "title" in info: # this can only have been set via menuItem ^^^
                deprecated("[config file] coloring '{}': 'legendTitle' has been replaced with 'title' & is unused since 'menuItem' is present.".format(key))
            else:
                deprecated("[config file] coloring '{}': 'legendTitle' has been replaced with 'title'. Using 'legendTitle' as 'title'.".format(key))
                info["title"] = info["legendTitle"]
            del info["legendTitle"]
        if "key" in info:
            del info["key"]
        if info.get("type") == "discrete":
            deprecated("[config file] coloring '{}': type 'discrete' is no longer valid. Please use 'ordinal', 'categorical' or 'boolean'. "
                "This has been automatically changed to 'categorical'.".format(key))
            info["type"] = "categorical"
        colorings.append({"key": key, **info})
    return colorings

def _rename_deprecated_colorings(colorings: list):
    for old_key, new_key in DEPRECATED_KEYS.items():
        new_coloring  = next(iter([{'index': idx, 'coloring': c} for [idx, c] in enumerate(colorings) if c['key']==new_key] or [None]))
        old_coloring  = next(iter([{'index': idx, 'coloring': c} for [idx, c] in enumerate(colorings) if c['key']==old_key] or [None]))
        if old_coloring:
            if new_coloring:
                deprecated(f"[config file] ignoring deprecated {old_key!r} because a coloring for {new_key!r} already exists")
                colorings = [c for [idx, c] in enumerate(colorings) if idx!=old_coloring['index']]
            else:
                deprecated(f"[config file] renaming coloring {old_key!r} to {new_key!r}")
                colorings[old_coloring['index']]['key'] = new_key
    return colorings

def _rename_display_keys(display: dict) -> dict:
    defaults = {**display}
    v1_v2_keys = [
        # [ <old key>, <new key> ]
        ["geoResolution", "geo_resolution"],
        ["colorBy", "color_by"],
        ["distanceMeasure", "distance_measure"],
        ["mapTriplicate", "map_triplicate"]
    ]
    for [v1_key, v2_key] in [x for x in v1_v2_keys if x[0] in defaults]:
        deprecated("[config file] '{}' has been replaced with '{}'".format(v1_key, v2_key))
        defaults[v2_key] = defaults[v1_key]
        del defaults[v1_key]
    return defaults

def _rename_deprecated_filters(filters: list) -> list:
    return [update_deprecated_names(f) for f in filters]

def remove_unused_metadata_columns(columns: list) -> list:
    # 1. Remove any occurrences of 'author' and 'num_date' as it's a no-op - author
    #    and numerical date information is always exported on nodes if it's available.
    # 2. Remove "authors" as that was historically corrected to "author" during parsing,
    #    and similarly for 'numdate'
    return [n for n in columns if n not in ['author', 'authors', 'num_date', 'numdate']]

DEPRECATIONS = [
    {"old_name": 'vaccine_choices', "new_name": None}, # removed from config
    {"old_name": "updated", "new_name": None}, # removed from config
    {"old_name": "defaults", "new_name": "display_defaults", "modify": _rename_display_keys},
    {"old_name": "maintainer", "new_name": "maintainers", "modify": lambda m: [{"name": m[0], "url": m[1]}]},
    {"old_name": "geo", "new_name": "geo_resolutions", "modify": lambda values: [{"key": v} for v in values]},
    {"old_name": "color_options", "new_name": "colorings", "modify": _parse_color_options},
    {"old_name": "colorings", "new_name": "colorings", "modify": _rename_deprecated_colorings},
    {"old_name": "filters", "new_name": "filters", "modify": _rename_deprecated_filters},
    {"old_name": "metadata_columns", "new_name": "metadata_columns", "modify": remove_unused_metadata_columns},
]


def read_single_auspice_config(fname: str, validation_mode: ValidationMode) -> dict[str, Any]:
    config = read_json(fname)

    if validation_mode is not ValidationMode.SKIP:
        try:
            print(f"Validating config file {fname!r} against the JSON schema")
            validate_auspice_config_v2(fname)
        except ValidateError:
            print(f"Validation of {fname!r} failed. Please check the formatting of this file & refer to the augur documentation for further help. ")
            validation_failure(validation_mode)

    for data in DEPRECATIONS:
        _replace_deprecated(config, **data) # type: ignore[arg-type]

    return config

def merge_configs(configs: list[dict[str, Any]], validation_mode: ValidationMode, output_fname: Union[str,None]=None) -> dict[str, Any]:
    merged = configs[0]
    for overlay in configs[1:]:
        # We could leverage the config schemas here in the future if desired
        merged = _merge_scalar(merged, overlay, 'title')
        merged = _merge_lists(merged, overlay, 'colorings', lambda x,y: x.get('key')==y.get('key'))
        merged = _merge_lists(merged, overlay, 'geo_resolutions', lambda x,y: _geo_resolution_id(x)==_geo_resolution_id(y))
        merged = _merge_lists(merged, overlay, 'maintainers', lambda x,y: x.get('name')==y.get('name'))
        merged = _merge_scalar(merged, overlay, 'build_url')
        merged = _merge_scalar(merged, overlay, 'build_avatar')
        merged = _merge_lists(merged, overlay, 'filters', lambda x,y: x==y)
        merged = _merge_dicts(merged, overlay, 'display_defaults')
        merged = _merge_lists(merged, overlay, 'panels', lambda x,y: x==y)
        merged = _merge_lists(merged, overlay, 'data_provenance', lambda x,y: x.get('name')==y.get('name'))
        merged = _merge_lists(merged, overlay, 'metadata_columns', lambda x,y: x==y)

        # extensions have any type (as per the schema), but in practice I've only seen dicts used.
        # We merge these by taking the first encountered extension block, and if it's a dict
        # then any future extension blocks which are also dicts are merged in
        if 'extensions' not in merged and 'extensions' in overlay:
            merged['extensions'] = overlay['extensions']
        elif isinstance(merged.get('extensions', None), dict) and isinstance(overlay.get('extensions', None), dict):
            merged = _merge_dicts(merged, overlay, 'extensions')

    # Write the merged config before validation (if requested)
    if output_fname:
        print(f"Writing merged auspice config JSON to {output_fname!r}")
        # don't use our util function `write_json` as we don't want any modification of values
        with open_file(output_fname, 'w', encoding='utf-8') as handle:
            json.dump(merged, handle, indent=2)

    if validation_mode is not ValidationMode.SKIP:
        try:
            print(f"Validating merged config file against the JSON schema")
            validate_auspice_config_v2(merged)
        except ValidateError:
            if output_fname:
                print(f"Validation of merged config file failed. The merged JSON has been written to {output_fname!r} for debugging purposes", file=sys.stderr)
            else:
                print(f"Validation of merged config file failed. Here is the contents of the merged config JSON for debugging purposes:", file=sys.stderr)
                json.dump(merged, sys.stderr, indent=2)
            validation_failure(validation_mode)

    return merged


def read_auspice_configs(*fnames: str, validation_mode: ValidationMode, output_fname: Union[str,None]) -> dict[str, Any]:
    """
    Parses zero, one or multiple files as auspice config JSONs and returns a (merged) config dict.
    """
    configs = [read_single_auspice_config(fname, validation_mode) for fname in fnames]

    if len(configs)==0: # no config JSONs provided
        return {}
    elif len(configs)==1: # no merging necessary
        return configs[0]
    return merge_configs(configs, validation_mode, output_fname)
