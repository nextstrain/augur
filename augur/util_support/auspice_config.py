import sys
import os
import json
from augur.io.file import open_file
from collections import defaultdict
from typing import Union
from collections.abc import Callable

def _read_json(fname):
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

def deprecated(message):
    print("[deprecation TODO XXX]", message)
    # warn(message, DeprecationWarning, stacklevel=2)
    # global deprecationWarningsEmitted
    # deprecationWarningsEmitted=True

def _replace_deprecated(config: dict, old_name: str, new_name: str, modify: Union[None, Callable]):
    if new_name in config and old_name in config:
        deprecated(f"[config file] top level keys {new_name!r} and {old_name!r} were both present, " \
                   "however the latter is a deprecated version of the former and will be ignored.")
        del config[old_name]
    elif old_name in config:
        deprecated(f"[config file] top level key {old_name!r} is deprecated. Renaming to {new_name!r}" +
                   (" (with modifications)" if modify else ""))
        config[new_name] = modify(config[old_name]) if modify else config[old_name]
        del config[old_name]

def _remove_deprecated(config: dict, key: str):
    if key in config:
        deprecated(f"[config file] top level key {key!r} is no longer used.")
        del config[key]

def _parse_color_options(options: dict) -> list[dict]:
    # parse v1-style 'color_options' & convert to v2-style 'colorings'
    colorings = []
    for key, info in options.items():
        # note: if both legentTitle & menuItem are present then we use the latter. See https://github.com/nextstrain/auspice/issues/730
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
            deprecated("[config file] coloring '{}': type 'discrete' is no longer valid. Please use either 'ordinal', 'categorical' or 'boolean'. "
                "This has been automatically changed to 'categorical'.".format(key))
            info["type"] = "categorical"
        colorings.append({"key": key, **info})
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


TOP_LEVEL_DEPRECATED_KEYS = [
    # [ <old key>, <new key>, <modify func>] #
    ["defaults", "display_defaults", _rename_display_keys],
    ["maintainer", "maintainers", lambda m: [{"name": m[0], "url": m[1]}]],
    ["geo", "geo_resolutions", lambda values: [{"key": v} for v in values]],
    ["color_options", "colorings", _parse_color_options]
]

TOP_LEVEL_UNUSED_KEYS = ['vaccine_choices', 'updated']

def read_auspice_config(fname):

    config = _read_json(fname)

    for keys in TOP_LEVEL_DEPRECATED_KEYS:
        _replace_deprecated(config, *keys)
    for key in TOP_LEVEL_UNUSED_KEYS:
        _remove_deprecated(config, key)

    return config
