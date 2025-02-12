import sys
import os
import json
from ..io.file import open_file
from ..errors import AugurError
from .warnings import deprecated
from collections import defaultdict
from typing import Union, Any
from collections.abc import Callable


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

def _replace_deprecated(config: dict[str,Any], old_name: str, new_name: Union[None,str], modify: Union[None, Callable]=None):
    if old_name not in config:
        return # NO-OP
    elif new_name is None:
        deprecated(f"[config file] key {old_name!r} is no longer used and has been dropped from your config")
        del config[old_name]
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


DEPRECATIONS = [
    {"old_name": 'vaccine_choices', "new_name": None}, # removed from config
    {"old_name": "updated", "new_name": None}, # removed from config
    {"old_name": "defaults", "new_name": "display_defaults", "modify": None},
    {"old_name": "maintainer", "new_name": "maintainers", "modify": lambda m: [{"name": m[0], "url": m[1]}]},
    {"old_name": "geo", "new_name": "geo_resolutions", "modify": lambda values: [{"key": v} for v in values]},
    {"old_name": "color_options", "new_name": "colorings", "modify": _parse_color_options},
]


def read_auspice_config(fname):

    config = read_json(fname)

    for data in DEPRECATIONS:
        _replace_deprecated(config, **data)

    return config
