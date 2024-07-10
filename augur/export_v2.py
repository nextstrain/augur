"""
Export version 2 JSON schema for visualization with Auspice
"""
import os
from pathlib import Path
import sys
import time
from collections import defaultdict, deque, OrderedDict
import warnings
import numbers
import math
import re
from Bio import Phylo
from typing import Dict, Union, TypedDict, Any, Tuple

from .argparse_ import ExtendOverwriteDefault
from .errors import AugurError
from .io.file import open_file
from .io.metadata import DEFAULT_DELIMITERS, DEFAULT_ID_COLUMNS, InvalidDelimiter, read_metadata
from .types import ValidationMode
from .utils import read_node_data, write_json, json_size, read_config, read_lat_longs, read_colors
from .validate import export_v2 as validate_v2, auspice_config_v2 as validate_auspice_config_v2, ValidateError


MINIFY_THRESHOLD_MB = 5


# Set up warnings & exceptions
warn = warnings.warn
deprecationWarningsEmitted = False

def deprecated(message):
    warn(message, DeprecationWarning, stacklevel=2)
    global deprecationWarningsEmitted
    deprecationWarningsEmitted=True

def warning(message):
    warn(message, UserWarning, stacklevel=2)

def fatal(message):
    print("FATAL ERROR: {}".format(message))
    sys.exit(2)

def configure_warnings():
    # we must only set these when someone runs `augur export v2` (i.e. run_v2() is called)
    # else they will apply to _all_ augur commands due to the way all commands are pulled
    # in by the augur runner (augur/__init__.py)
    def customformatwarning(message, category, filename, lineno, line=None):
        if category.__name__ == "UserWarning":
            return "WARNING: {}\n\n".format(message)
        if category.__name__ == "DeprecationWarning":
            return "DEPRECATED: {}\n\n".format(message)
        return "{}\n".format(message)

    warnings.formatwarning = customformatwarning
    warnings.simplefilter("default") # show DeprecationWarnings by default

class InvalidOption(Exception):
    pass

class CustomOrderedDict(OrderedDict):
    """
    Similar to OrderedDict but will convert dictionaries (and dictionaries of dictionaries)
    into (nested) CustomOrderedDicts.
    Encountered lists of dicts will be converted to lists of CustomOrderedDict but we will not
    recursively explore nested lists.
    Tuples and other iterators are not explored.
    """
    def __init__(self, *args):
        super().__init__(*args)
        for key in self:
            if isinstance(self[key], dict) and not isinstance(self[key], OrderedDict):
                self[key] = CustomOrderedDict(self[key])
            elif isinstance(self[key], list):
                self[key] = [
                    (CustomOrderedDict(el) if (isinstance(el, dict) and not isinstance(el, OrderedDict)) else el)
                    for el in self[key]
                ]
    def set_order(self, *order):
        """
        changes the order of keys to match those specified in `order` as much
        as possible. Missing keys are ignored. Extra keys will come after those
        specified in `order`.
        """
        for key in reversed(order):
            self.move_to_end_if_present(key, last=False)
    def move_to_end_if_present(self, key, **kwargs):
        try:
            self.move_to_end(key, **kwargs)
        except KeyError:
            pass


def orderKeys(data):
    """
    converts the data dict (where keys are inherently unordered) into an
    OrderedDict where keys are nicely ordered for human eyes to scan the
    data when written to JSON. The ordering (mostly) mirrors the schema.
    """
    od = CustomOrderedDict(data)
    od.set_order("version", "meta", "tree")
    if "meta" in od:
        od["meta"].set_order("title", "updated", "build_url", "data_provenance", "maintainers")
        for coloring in od['meta'].get('colorings', []):
            coloring.set_order("key", "title", "type", "scale", "legend")
    def order_nodes(node):
        """recursive function to order nodes in a (sub)tree"""
        node.set_order("name", "node_attrs", "branch_attrs")
        # children often a _large_ object and it improves readability if this comes last in the node
        node.move_to_end_if_present("children")
        if "node_attrs" in node:
            node["node_attrs"].set_order("div", "num_date")
        for child in node.get("children", []):
            order_nodes(child)
    if isinstance(od.get("tree"), list):
        for subtree in od['tree']:
            order_nodes(subtree)
    elif isinstance(od.get("tree"), dict):
        order_nodes(od['tree'])
    return od


def node_div(T, node_attrs):
    """
    Scans the provided tree & metadata to see if divergence is defined, and if so returns
    a function which gets it from individual nodes. Divergence may be defined via a number
    of sources, and we pick them in the following order:
    * metadata.mutation_length (typically via `augur refine`)
    * metadata.branch_length (typically via `augur refine`)
    * Branch lengths encoded in the Newick tree

    Returns either:
    * function with arguments: (node, metadata_for_node) which returns the node divergence
    * None (indicates that divergence is not available for this dataset)
    """
    if all(('mutation_length' in node_attrs[n.name] for n in T.root.find_clades())):
        return lambda node, metadata: metadata['mutation_length']
    if all(('branch_length' in node_attrs[n.name] for n in T.root.find_clades())):
        return lambda node, metadata: metadata['branch_length']
    if T.root.branch_length is not None:
        return lambda node, metadata: node.branch_length
    return None

def convert_tree_to_json_structure(node, metadata, get_div, div=0):
    """
    converts the Biopython tree structure to a dictionary that can
    be written to file as a json. This is called recursively.

    Parameters
    ----------
    node : Bio.Phylo.Newick.Clade
    metadata : dict
        Per-node metadata, with keys matching `node.name`
    get_div :
        (None or function)
        Function returns divergence for this node. Arguments: (node, metadata_for_node)
        If None then divergence is not defined for this dataset and so 'div' is not set on returned nodes.
    div : int
        cumulative divergence leading to the current node (root = 0)

    Returns
    -------
    dict:
        See schema-export-v2.json#/$defs/tree for full details.
        Node names are always set, and divergence is set if applicable
    """
    node_struct = {'name': node.name, 'node_attrs': {}, 'branch_attrs': {}}

    if get_div is not None: # Store the (cumulative) observed divergence prior to this node
        node_struct["node_attrs"]["div"] = format_number(div)

    if node.clades:
        node_struct["children"] = []
        for child in node.clades:
            cdiv = div
            if get_div:
                cdiv += get_div(child, metadata[child.name])
            node_struct["children"].append(convert_tree_to_json_structure(child, metadata, get_div, div=cdiv))

    return node_struct

def are_mutations_defined(branch_attrs):
    for branch in branch_attrs.values():
        if branch.get("mutations"):
            return True
    return False

def is_node_attr_defined(node_attrs, attr_name):
    for node, data in node_attrs.items():
        if data.get(attr_name):
            return True
    return False


def update_deprecated_names(name):
    # correct deprecated keys
    change = {
        "authors": "author",
        "numdate": "num_date"
    }
    return change.get(name, name)

def get_values_across_nodes(node_attrs, key):
    vals = set()
    for data in node_attrs.values():
        if is_valid(data.get(key)):
            vals.add(data.get(key))
    return vals

def is_valid(value):
    invalid = ["undefined", "unknown", "?", "nan", "na", "n/a", 'none', '', 'not known']
    return str(value).strip('\"').strip("'").strip().lower() not in invalid

def get_config_colorings_as_dict(config):
    # extract the colorings information from the config JSON, if provided
    # and return as a dict of `key` -> obj where `obj` is in v2-schema format
    config_colorings = {}
    if config.get("colorings"):
        config_colorings = {v.get("key"):v for v in config["colorings"]}
    elif config.get("color_options"):
        deprecated("[config file] 'color_options' has been replaced with 'colorings' & the structure has changed.")
        # parse v1-style colorings & convert to v2-style
        for key, info in config.get("color_options").items():
            # note: if both legentTitle & menuItem are present then we use the latter. See https://github.com/nextstrain/auspice/issues/730
            if "menuItem" in info:
                deprecated("[config file] coloring '{}': 'menuItem' has been replaced with 'title'. Using 'menuItem' as 'title'.".format(key))
                info["title"] = info["menuItem"]
                del info["menuItem"]
            if "legendTitle" in info:
                if info["title"]: # this can only have been set via menuItem ^^^
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
            config_colorings[key] = info
    return config_colorings


def set_colorings(data_json, config, command_line_colorings, metadata_names, node_data_colorings, provided_colors, node_attrs, branch_attrs):

    def _get_type(key, trait_values):
        # for some keys we know what the type must be
        known_types = {
            "clade_membership": "categorical",
            "gt": "categorical",
            "author": "categorical",
            "num_date": "continuous"
        }
        if key in known_types:
            return known_types[key]

        if config.get(key, {}).get("type"):
            t = config.get(key).get("type")
            allowedTypes = ["continuous", "temporal", "ordinal", "categorical", "boolean"]
            if t not in allowedTypes:
                warn("[config file] In trait {}, coloring type '{}' is not valid. Please choose from: '{}'. This trait has been excluded!".format(key, t, ", ".join(allowedTypes)))
                raise InvalidOption()
            return t
        # no type supplied => try to guess
        if all([all([str(x).lower() in ["false", "true", "1.0", "0.0", "1", "0", "yes", "no"] for x in trait_values])]):
            t = "boolean"
        elif all([ isinstance(n, float) if isinstance(n, float) else isinstance(n, int) for n in trait_values ]):
            t = "continuous"
        else:
            t = "categorical"
        #Don't warn if command-line - no way to specify
        #treat country and region differently
        if config:
            warn("[config file] Trait '{}' is missing type information. We've guessed '{}'.".format(key, t))
        elif key != "country" and key != "region":
            print("Trait '{}' was guessed as being type '{}'. Use a 'config' file if you'd like to set this yourself.".format(key, t))
        return t

    def _get_title(key):
        # preferentially get the title from the config if set (walrus operator would help here)
        config_title = config.get(key, {}).get("title")
        if config_title:
            return config_title

        # hardcoded fallbacks:
        if key == "clade_membership":
            return "Clade"
        if key == "gt":
            return "Genotype"
        if key == "author":
            return "Authors"
        if key == 'num_date':
            return 'Sampling date'
        # fallthrough
        return key

    def _add_color_scale(coloring):
        ## consider various sources to find any user-provided scale information
        key = coloring["key"]
        scale_type = coloring["type"]
        if scale_type == "continuous":
            ## continuous scale information can only come from an auspice config JSON
            if config.get(key, {}).get("scale"):
                # enforce numeric values (we can't use the schema for this)
                provided_scale = [s for s in config[key]['scale'] if isinstance(s[0], numbers.Number)]
                if len(provided_scale)<2:
                    warn(f"The scale provided for {key} had fewer than 2 (valid numeric) entries. Skipping.")
                    return coloring
                coloring["scale"] = provided_scale
                return coloring
        elif config.get(key, {}).get("scale"):
            # If the auspice config JSON (`config`) explicitly defines a scale for this coloring
            # then we use this instead of any colors provided via a TSV file (`provided_colors`)
            values_in_tree = get_values_across_nodes(node_attrs, key)
            scale = []
            provided_values_unseen_in_tree = []
            for info in config[key]['scale']:
                if info[0] in values_in_tree:
                    scale.append(info)
                else:
                    provided_values_unseen_in_tree.append(info[0])
            if len(scale):
                coloring["scale"] = scale
                if len(provided_values_unseen_in_tree):
                    warn(f"The configuration JSON specifies colors for \"{key}\" which aren't found in the tree:\n\t{', '.join(provided_values_unseen_in_tree)}.")
                return coloring
            warn(f"The configuration JSON specifies a color scale for {key} however none of the values in the tree are in this scale! Auspice will generate its own color scale for this trait.")
        elif key.lower() in provided_colors:
            # `provided_colors` typically originates from a colors.tsv file
            scale = []
            trait_values = {str(val).lower(): val for val in get_values_across_nodes(node_attrs, key)}
            trait_values_unseen = {k for k in trait_values}
            for provided_key, provided_color in provided_colors[key.lower()]:
                if provided_key.lower() in trait_values:
                    scale.append([trait_values[provided_key.lower()], provided_color])
                    trait_values_unseen.discard(provided_key.lower())
            if len(scale):
                coloring["scale"] = scale
                if len(trait_values_unseen):
                    warn(f"These values for trait {key} were not specified in the colors file you provided:\n\t{', '.join(trait_values_unseen)}.\n\tAuspice will create colors for them.")
                return coloring
            warn(f"You've supplied a colors file with information for {key} but none of the values found on the tree had associated colors. Auspice will generate its own color scale for this trait.")
        # Fallthrough (no scale information provided means that auspice will generate its own scale)
        return coloring

    def _add_legend(coloring):
        """
        If there is a config-defined legend entry for this coloring, then add it to
        the coloring object under the key "legend"

        Enclosing scope variables used
        ------------------------------
        config : object

        Parameters
        ----------
        coloring : object

        Returns
        -------
        object :
            returns the single input parameter, potentially modified in place
        """
        key = coloring["key"]
        if config.get(key, {}).get("legend"):
            # The structure of the legend in the auspice config file is the same as the exported
            # auspice dataset JSON, and the schema validates this for us so we don't have to check here.
            # Note that if the data is inconsistent (e.g. overlapping bounds) then auspice will
            # discard them (and print a console warning)
            coloring['legend'] = config[key]['legend']
        return coloring

    def _add_title_and_type(coloring):
        key = coloring["key"]
        trait_values = get_values_across_nodes(node_attrs, key) # e.g. list of countries, regions etc
        try:
            coloring["title"] = _get_title(key)
            coloring["type"] = _get_type(key, trait_values)
        except InvalidOption:
            return False # a warning message will have been printed before `InvalidOption` is raised
        return coloring

    def _add_coloring(colorings, key):
        # handle deprecations
        if key == "authors":
            deprecated("[colorings] The 'authors' key is now called 'author'")
            key = "author"
        # check if the key has already been added by another part of the color-creating logic
        if key not in {x['key'] for x in colorings}:
            colorings.append({"key": key})

    def _is_valid(coloring):
        key = coloring["key"]
        trait_values = get_values_across_nodes(node_attrs, key) # e.g. list of countries, regions etc
        if key == "gt" and not are_mutations_defined(branch_attrs):
            warn("[colorings] You asked for mutations (\"gt\"), but none are defined on the tree. They cannot be used as a coloring.")
            return False
        if key != "gt" and not trait_values:
            warn(f"Requested color-by field {key!r} does not exist and will not be used as a coloring or exported.")
            return False
        return True

    def _get_colorings():
        # which colorings are intended for this dataset?
        # returns an array of dicts, where the order determines the order in auspice's dropdowns
        # note that invalid options will be pruned out later
        # it is here that we deal with the interplay between node-data "traits", command line colorings &
        # config provided options
        auto_colorings = [name for name in node_data_colorings
                          if node_data_prop_is_normal_trait(name) and name not in metadata_names]

        colorings = []
        # If we have command line colorings, it seems we (a) ignore any provided in the config file
        # and (b) start with the node-data "traits". (Note that in a later function, the title and/or
        # type will be accessed from the config file if available)
        if command_line_colorings:
            # start with auto_colorings (already validated to be included)
            for x in auto_colorings:
                _add_coloring(colorings, x)
            # then add in command line colorings
            for x in command_line_colorings:
                _add_coloring(colorings, x)
        else:
            # if we have a config file, start with these (extra info, such as title&type, is added in later)
            if config:
                for x in config.keys():
                    _add_coloring(colorings, x)
            # then add in any auto-colorings already validated to include
            for x in auto_colorings:
                _add_coloring(colorings, x)

        explicitly_defined_colorings = [x["key"] for x in colorings]
        # add in genotype as a special case if (a) not already set and (b) the data supports it
        if "gt" not in explicitly_defined_colorings and are_mutations_defined(branch_attrs):
            colorings.insert(0,{'key':'gt'})
        if "num_date" not in explicitly_defined_colorings and is_node_attr_defined(node_attrs, "num_date"):
            colorings.insert(0,{'key':'num_date'})
        if "clade_membership" not in explicitly_defined_colorings and is_node_attr_defined(node_attrs, "clade_membership"):
            colorings.insert(0,{'key':'clade_membership'})

        return colorings


    # construct colorings from cmd line args, data, config file etc
    colorings = _get_colorings()
    # ensure the data supports each coloring & ignore if not
    colorings = [c for c in colorings if _is_valid(c)]
    # add the title / type from config or via predefined logic rules. Note this can return False on an error.
    colorings = [x for x in [_add_title_and_type(coloring) for coloring in colorings] if x]
    # for each coloring, if colors have been provided, save them as a "scale"
    colorings = [_add_color_scale(coloring) for coloring in colorings]
    # for each coloring, pass-through the `legend` data from auspice_config.json, if provided
    colorings = [_add_legend(coloring) for coloring in colorings]
    # save them to the data json to be exported
    data_json['meta']["colorings"] = colorings


def set_geo_resolutions(data_json, config, command_line_traits, lat_long_mapping, node_attrs):
    """
    appropriately combine provided geo resolutions from command line & config files
    and associate with lat/longs.
    """
    geo_resolutions = []

    def _transfer_geo_data(node):
        for g in geo_resolutions:
            if g['key'] in node_attrs[node["name"]] and g['key'] not in node['node_attrs'] \
                and is_valid(node_attrs[node["name"]][g['key']]): # don't add if not valid!
                node['node_attrs'][g['key']] = {"value":node_attrs[node["name"]][g['key']]}

        if "children" in node:
            for c in node["children"]:
                _transfer_geo_data(c)

    # step 1: get a list of resolutions
    if command_line_traits:
        # straight overwrite -- not an extension of those which may be provided in the config
        traits = [{"key": x} for x in command_line_traits]
    elif config.get("geo_resolutions"):
        config_geo = config.get("geo_resolutions")
        # If is set up correctly as dict with "key", take it as-is
        if all( (isinstance(entry, dict) and "key" in entry.keys()) for entry in config_geo):
            traits = config.get("geo_resolutions")
        # If is a list of strings, or mix of strings and dicts with "key", we can do this!
        elif all( (isinstance(entry, dict) and "key" in entry.keys()) or isinstance(entry, str) for entry in config_geo):
            traits = [entry if isinstance(entry, dict) else {"key": entry} for entry in config_geo]
        else:
            print("WARNING: [config file] 'geo_resolutions' is not in an acceptible format. The field is now list of strings, or list of dicts each with format {\"key\":\"country\"}")
            print("\t It is being ignored - this run will have no geo_resolutions.\n")
            return False
    elif config.get("geo"):
        traits = [{"key": x} for x in config.get("geo")]
        deprecated("[config file] 'geo' has been replaced with 'geo_resolutions'. The field is now list of strings, or list of dicts each with format {\"key\":\"country\"}")
    else:
        return False

    # step 2: for each resolution, create the map of deme name -> lat/long
    for trait_info in traits:
        deme_to_lat_longs = {}
        trait_values = get_values_across_nodes(node_attrs, trait_info["key"]) # e.g. list of countries, regions etc

        for deme in trait_values: # note: deme may be numeric, or string
            try:
                deme_search_value = deme.lower()
            except AttributeError:
                deme_search_value = str(deme)
            try:
                deme_to_lat_longs[deme] = lat_long_mapping[(trait_info["key"].lower(), deme_search_value)]
            except KeyError:
                warn("{}->{} did not have an associated lat/long value (matching performed in lower case). Auspice won't be able to display this location.".format(trait_info["key"], deme))

        if deme_to_lat_longs:
            data = {"key": trait_info["key"], "demes": deme_to_lat_longs}
            if "title" in trait_info:
                data["title"] = trait_info["title"]
            geo_resolutions.append(data)
        else:
            warn("Geo resolution \"{}\" had no demes with supplied lat/longs and will be excluded from the exported \"geo_resolutions\".".format(trait_info["key"]))

    #
    _transfer_geo_data(data_json['tree'])


    if geo_resolutions:
        data_json['meta']["geo_resolutions"] = geo_resolutions


def set_annotations(data_json, node_data):
    if "annotations" in node_data:
        data_json['meta']["genome_annotations"] = node_data["annotations"]

def set_filters(data_json, config):
    # NB set_colorings() must have been run as we access those results
    if config.get('filters') == []:
        # an empty config section indicates no filters are to be exported
        data_json['meta']['filters'] = []
    elif config.get('filters'):
        # set filters as long as they are not continuous
        data_json['meta']['filters'] = config['filters']
    else:
        # if not specified in the config, include all boolean and categorical colorbys
        potentials = {coloring["key"] for coloring in data_json['meta']["colorings"]
                      if coloring["type"] != "continuous" and coloring["key"] != 'gt'}
        data_json['meta']['filters'] = list(potentials)

def validate_data_json(filename, validation_mode=ValidationMode.ERROR):
    if validation_mode is ValidationMode.SKIP:
        print(f"Skipping validation of produced JSON due to --validation-mode={validation_mode.value} or --skip-validation.")
        return

    print("Validating produced JSON")
    try:
        validate_v2(main_json=filename)
    except ValidateError as e:
        print(e)
        print("\n------------------------")
        print("Validation of {} failed. Please check this in a local instance of `auspice`, as it is not expected to display correctly. ".format(filename))
        print("------------------------")
        validation_failure(validation_mode)

def validation_failure(mode: ValidationMode):
    if mode is ValidationMode.ERROR:
        sys.exit(2)
    elif mode is ValidationMode.WARN:
        print(f"Continuing due to --validation-mode={mode.value} even though there were validation errors.")
    elif mode is ValidationMode.SKIP:
        # Shouldn't be doing validation under skip, but if we're called anyway just do nothing.
        return
    else:
        raise ValueError(f"unknown validation mode: {mode!r}")


def set_panels(data_json, config, cmd_line_panels):

    # config set panels overrides cmd-line provided panels
    panels = config["panels"] if config.get("panels") else cmd_line_panels

    if not panels:
        panels = ["tree", "map", "entropy"]

    if data_json["meta"].get("genome_annotations", None):
        annotations = data_json["meta"]["genome_annotations"].keys()
    else:
        annotations = []

    if "entropy" in panels and len(annotations) == 0:
        panels.remove("entropy")

    # for map to be displayed, we need to have valid geo resolutions
    if "map" in panels:
        if "geo_resolutions" not in data_json["meta"] or not data_json["meta"]["geo_resolutions"]:
            panels.remove("map")

    data_json['meta']["panels"] = panels


def set_data_provenance(data_json, config):
    """Set the data provenance from the given config file to the given data JSON.

    Parameters
    ----------
    data_json : dict
        auspice JSON to be updated
    config : dict
        config JSON with an expected ``data_provenance`` key

    Examples
    --------
    >>> config = {"data_provenance": [{"name": "GISAID"}, {"name": "INSDC"}]}
    >>> data_json = {"meta": {}}
    >>> set_data_provenance(data_json, config)
    >>> data_json["meta"]["data_provenance"][0]["name"]
    'GISAID'

    """
    if "data_provenance" in config:
        data_json["meta"]["data_provenance"] = config["data_provenance"]


def counter_to_disambiguation_suffix(count):
    """Given a numeric count of author papers, return a distinct alphabetical
    disambiguation suffix.

    Examples
    --------
    >>> counter_to_disambiguation_suffix(0)
    'A'
    >>> counter_to_disambiguation_suffix(25)
    'Z'
    >>> counter_to_disambiguation_suffix(26)
    'AA'
    >>> counter_to_disambiguation_suffix(51)
    'AZ'
    >>> counter_to_disambiguation_suffix(52)
    'BA'
    """
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    base = len(letters)
    suffix = deque()

    # Find the appropriate combination of letters for the given count. This
    # closely resembles the steps required to calculate the base 26 value of the
    # given base 10 number.
    while True:
        quotient = count // base
        remainder = count % base

        # Collect remainders from right to left. Letters are zero-indexed such
        # that a count of 0 returns an "A".
        suffix.appendleft(letters[remainder])

        # Stop when we've accounted for all possible quotient and remainder
        # values.
        if quotient == 0:
            break

        # Convert counts to zero-indexed values such that the next place value
        # starts with the letter "A" instead of the letter "B".
        count = quotient - 1

    return "".join(suffix)

def create_author_data(node_attrs):
    """Gather the authors which appear in the metadata and create the author
    info structure with unique keys
    """

    def node_to_author_tuple(data):
        # make a unique list of citations to disambiguate
        # Hadfield et al A, Hadfield et al B, etc...
        return (
            data.get("author", "unknown"),
            data.get("title", "unknown"),
            data.get("journal", "unknown")
        )

    author_to_unique_tuples = defaultdict(list)
    node_author_info = {}

    # first pass:
    for node_name, node_info in node_attrs.items():
        author = node_info.get("author")
        if not author:
            author = node_info.get("authors")
        if not author:
            continue # internal node / terminal node without authors

        node_author_info[node_name] = {"author": author}
        if "title" in node_info:
            title = node_info["title"].strip()
            if is_valid(title):
                node_author_info[node_name]["title"] = title
        if "journal" in node_info:
            journal = node_info["journal"].strip()
            if is_valid(journal):
                node_author_info[node_name]["journal"] = journal
        if "paper_url" in node_info:
            paper_url = node_info["paper_url"].strip()
            if is_valid(paper_url):
                node_author_info[node_name]["paper_url"] = paper_url

        author_tuple = node_to_author_tuple(node_author_info[node_name])

        # relies on defaultdict to initialize empty list
        if author_tuple not in author_to_unique_tuples[author]:
            author_to_unique_tuples[author].append(author_tuple)

    # second pass is necessary because we don't know if we need A, B, C
    # without a complete first pass
    for node_name, node_info in node_attrs.items():
        if node_name not in node_author_info:
            continue # internal node / terminal node without authors

        author_tuple = node_to_author_tuple(node_author_info[node_name])
        author = node_author_info[node_name]["author"]
        if len(author_to_unique_tuples[author]) > 1:
            index = author_to_unique_tuples[author].index(author_tuple)
            disambiguation_suffix = counter_to_disambiguation_suffix(index)
            node_author_info[node_name]["value"] = f"{author} {disambiguation_suffix}"
        else:
            node_author_info[node_name]["value"] = author

    return node_author_info

def set_branch_attrs_on_tree(data_json, branch_attrs):
    """
    Shifts the provided `branch_attrs` onto the (auspice) `data_json`.
    Currently all data is transferred, there is no way for (e.g.) the set of exported
    labels to be restricted by the user in a config.
    """
    def _recursively_set_data(node):
        if branch_attrs.get(node['name'], {}):
            node['branch_attrs'] = branch_attrs[node['name']]
        for child in node.get("children", []):
            _recursively_set_data(child)
    _recursively_set_data(data_json["tree"])

def is_numeric(n:Any) -> bool:
    # Typing a number is surprisingly hard in python, and `number.Number`
    # doesn't work nicely with type hints. See <https://stackoverflow.com/a/73186301>
    return isinstance(n, (int, float))

def format_number(n: Union[int, float]) -> Union[int, float]:
    if isinstance(n, int) or n==0:
        return n
    # We want to use three sig figs for the fractional part of the float, while
    # preserving all the data in the integral (integer) part. We leverage the
    # fact that python floats (incl scientific notation) storing the shortest
    # decimal string that’s guaranteed to round back to x. Note that this means we
    # drop trailing zeros, so it's not _quite_ sig figs.
    integral = int(abs(n))
    significand = math.floor(math.log10(integral))+1 if integral!=0 else 0
    return float(f"{n:.{significand+3}g}")


class ConfidenceNumeric(TypedDict):
    # the python type is a tuple, but when serialised to JSON this is an array
    confidence: Tuple[Union[int,float], Union[int,float]]

class ConfidenceCategorical(TypedDict):
    confidence: Dict[str,Union[int,float]]
    entropy: float

class EmptyDict(TypedDict):
    """Empty dict for typing."""

def attr_confidence(
    node_name: str,
    attrs: dict,
    key: str,
) -> Union[EmptyDict, ConfidenceNumeric, ConfidenceCategorical]:
    """
    Extracts and formats the confidence & entropy keys from the provided node-data attrs
    If there is no confidence-related information an empty dict is returned.
    If the information appears incorrect / incomplete, a warning is printed and an empty dict returned.
    """
    conf_key = f"{key}_confidence"
    conf = attrs.get(conf_key, None)
    if conf is None:
        return {}

    if isinstance(conf, list):
        if len(conf)!=2:
            warn(f"[confidence] node {node_name!r} specifies {conf_key!r} as a list of {len(conf)} values, not 2. Skipping confidence export.")
            return {}
        return {"confidence": (format_number(conf[0]), format_number(conf[1]))}

    if isinstance(conf, dict):
        entropy = attrs.get(f"{key}_entropy", None)
        if not entropy or not is_numeric(entropy):
            warn(f"[confidence] node {node_name!r} includes a mapping of confidence values but not an associated numeric entropy value. Skipping confidence export.")
            return {}
        if not all([is_numeric(v) for v in conf.values()]):
            warn(f"[confidence] node {node_name!r} includes a mapping of confidence values but they are not all numeric. Skipping confidence export.")
            return {}
        # While most of the time confidences come from `augur traits` which already sorts the values, we sort them (again) here
        # and only take confidence values over .1% and the top 4 elements.
        # To minimise the JSON size we only print values to 3 s.f. which is enough for Auspice
        conf = {
            key:format_number(conf[key]) for key in
            sorted(list(conf.keys()), key=lambda x: conf[x], reverse=True)
            if conf[key]>0.001
        }
        return {"confidence": conf, "entropy": format_number(entropy)}

    warn(f"[confidence] {key+'_confidence'!r} is of an unknown format. Skipping.")
    return {}



def set_node_attrs_on_tree(data_json, node_attrs, additional_metadata_columns):
    '''
    Assign desired colorings, metadata etc to the `node_attrs` of nodes in the tree

    Parameters
    ----------
    data_json : dict
    node_attrs: dict
        keys: strain names. values: dict with keys -> all available metadata (even "excluded" keys), values -> data (string / numeric / bool)
    additional_metadata_columns: list
        Requested additional metadata columns to export
    '''

    author_data = create_author_data(node_attrs)

    def _transfer_additional_metadata_columns(node, raw_data):
        for col in additional_metadata_columns:
            if is_valid(raw_data.get(col, None)):
                node["node_attrs"][col] = {"value": raw_data[col]}

    def _transfer_vaccine_info(node, raw_data):
        if raw_data.get("vaccine"):
            node["node_attrs"]['vaccine'] = raw_data['vaccine']

    def _transfer_hidden_flag(node, raw_data):
        hidden = raw_data.get("hidden", None)
        if hidden:
            if hidden in ["always", "divtree", "timetree"]:
                node["node_attrs"]["hidden"] = hidden
            elif hidden is True or str(hidden) == "1": # interpret this as hidden in both div + time tree
                node["node_attrs"]["hidden"] = "always"
            else:
                warn("Hidden node trait of {} is invalid. Ignoring.".format(hidden))

    def _transfer_num_date(node, raw_data):
        if raw_data.get("numdate", None) and not raw_data.get("num_date", None):
            raw_data["num_date"] = raw_data["numdate"]
            del raw_data["numdate"]
        if is_valid(raw_data.get("num_date", None)): # it's ok not to have temporal information
            node["node_attrs"]["num_date"] = {"value": format_number(raw_data["num_date"])}
            node["node_attrs"]["num_date"].update(attr_confidence(node["name"], raw_data, "num_date"))

    def _transfer_url_accession(node, raw_data):
        for prop in ["url", "accession"]:
            if is_valid(raw_data.get(prop, None)):
                node["node_attrs"][prop] = str(raw_data[prop])

    def _transfer_colorings_filters(node, raw_data):
        trait_keys = set() # order we add to the node_attrs is not important for auspice
        if "colorings" in data_json["meta"]:
            trait_keys = trait_keys.union([t["key"] for t in data_json["meta"]["colorings"]])
        if "filters" in data_json["meta"]:
            trait_keys = trait_keys.union(data_json["meta"]["filters"])
        exclude_list = ["gt", "num_date", "author"] # exclude special cases already taken care of
        trait_keys = trait_keys.difference(exclude_list)
        for key in trait_keys:
            value = raw_data.get(key, None)
            if is_valid(value):
                node["node_attrs"][key] = {"value": format_number(value) if is_numeric(value) else value}
                node["node_attrs"][key].update(attr_confidence(node["name"], raw_data, key))

    def _transfer_author_data(node):
        if node["name"] in author_data:
            node["node_attrs"]["author"] = author_data[node["name"]]

    def _recursively_set_data(node):
        # get all the available information for this particular node
        raw_data = node_attrs[node["name"]]
        # transfer requested metadata columns first so that the "special cases"
        # below can overwrite them as necessary
        _transfer_additional_metadata_columns(node, raw_data)
        # transfer "special cases"
        _transfer_vaccine_info(node, raw_data)
        _transfer_hidden_flag(node, raw_data)
        _transfer_num_date(node, raw_data)
        _transfer_url_accession(node, raw_data)
        _transfer_author_data(node)
        # transfer colorings & filters, including entropy & confidence if available
        _transfer_colorings_filters(node, raw_data)

        for child in node.get("children", []):
            _recursively_set_data(child)

    _recursively_set_data(data_json["tree"])

def node_data_prop_is_normal_trait(name):
    # those traits / keys / attrs which are not "special" and can be exported
    # as normal attributes on nodes
    excluded = [
        "authors",          # authors are set as a node property, not a trait property
        "author",           # see above
        "vaccine",          # vaccine info is stored as a "special" node prop
        'clade_membership', # explicitly set as a coloring if present
        'branch_length',
        'num_date',         # explicitly set as a coloring if present
        'raw_date',
        'numdate',
        'clock_length',
        'mutation_length',
        'date',
        'muts',
        'aa_muts',
        'sequence',
        'aa_sequences',
        'hidden',
        'dTiter',
        'dTiterSub'
    ]
    if name in excluded:
        return False

    if '_confidence' in name or '_entropy' in name:
        return False

    return True

validation_mode_help_message = """
    Control if optional validation checks are performed and what
    happens if they fail.

    'error' and 'warn' modes perform validation and emit messages about
    failed validation checks.  'error' mode causes a non-zero exit
    status if any validation checks failed, while 'warn' does not.

    'skip' mode performs no validation.

    Note that some validation checks are non-optional and as such are
    not affected by this setting.
"""


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("v2", help=__doc__)

    required = parser.add_argument_group(
        title="REQUIRED"
    )
    required.add_argument('--tree','-t', metavar="newick", required=True, help="Phylogenetic tree, usually output from `augur refine`")
    required.add_argument('--output', metavar="JSON", required=True, help="Output file (typically for visualisation in auspice)")

    config = parser.add_argument_group(
        title="DISPLAY CONFIGURATION",
        description="These control the display settings for auspice. \
            You can supply a config JSON (which has all available options) or command line arguments (which are more limited but great to get started). \
            Supplying both is fine too, command line args will overrule what is set in the config file!"
    )
    config.add_argument('--auspice-config', metavar="JSON", help="Auspice configuration file")
    config.add_argument('--title', type=str, metavar="title", help="Title to be displayed by auspice")
    config.add_argument('--maintainers', metavar="name", action="append", nargs='+', help="Analysis maintained by, in format 'Name <URL>' 'Name2 <URL>', ...")
    config.add_argument('--build-url', type=str, metavar="url", help="Build URL/repository to be displayed by Auspice")
    config.add_argument('--description', metavar="description.md", help="Markdown file with description of build and/or acknowledgements to be displayed by Auspice")
    config.add_argument('--geo-resolutions', metavar="trait", nargs='+', action='extend', help="Geographic traits to be displayed on map")
    config.add_argument('--color-by-metadata', metavar="trait", nargs='+', action='extend', help="Metadata columns to include as coloring options")
    config.add_argument('--metadata-columns', nargs="+", action="extend",
                                 help="Metadata columns to export in addition to columns provided by --color-by-metadata or colorings in the Auspice configuration file. " +
                                      "These columns will not be used as coloring options in Auspice but will be visible in the tree.")
    config.add_argument('--panels', metavar="panels", nargs='+', action='extend', choices=['tree', 'map', 'entropy', 'frequencies', 'measurements'], help="Restrict panel display in auspice. Options are %(choices)s. Ignore this option to display all available panels.")

    optional_inputs = parser.add_argument_group(
        title="OPTIONAL INPUT FILES"
    )
    optional_inputs.add_argument('--node-data', metavar="JSON", nargs='+', action="extend", help="JSON files containing metadata for nodes in the tree")
    optional_inputs.add_argument('--metadata', metavar="FILE", help="Additional metadata for strains in the tree")
    optional_inputs.add_argument('--metadata-delimiters', default=DEFAULT_DELIMITERS, nargs="+", action=ExtendOverwriteDefault,
                                 help="delimiters to accept when reading a metadata file. Only one delimiter will be inferred.")
    optional_inputs.add_argument('--metadata-id-columns', default=DEFAULT_ID_COLUMNS, nargs="+", action=ExtendOverwriteDefault,
                                 help="names of possible metadata columns containing identifier information, ordered by priority. Only one ID column will be inferred.")
    optional_inputs.add_argument('--colors', metavar="FILE", help="Custom color definitions, one per line in the format `TRAIT_TYPE\\tTRAIT_VALUE\\tHEX_CODE`")
    optional_inputs.add_argument('--lat-longs', metavar="TSV", help="Latitudes and longitudes for geography traits (overrides built in mappings)")

    minify_group = parser.add_argument_group(
            title="OPTIONAL MINIFY SETTINGS",
            description=f"""
                By default, output JSON files (both main and sidecar) are automatically minimized if
                the size of the un-minified main JSON file exceeds {MINIFY_THRESHOLD_MB} MB. Use
                these options to override that behavior.
                """
        ).add_mutually_exclusive_group()
    minify_group.add_argument('--minify-json', action="store_true", help="always export JSONs without indentation or line returns.")
    minify_group.add_argument('--no-minify-json', action="store_true", help="always export JSONs to be human readable.")

    root_sequence_group = parser.add_argument_group(
            title="OPTIONAL ROOT-SEQUENCE SETTINGS",
            description=f"""
                The root-sequences describe the sequences (nuc + aa) for the parent of the tree's root-node. They may represent a
                reference sequence or the inferred sequence at the root node, depending on how they were generated.
                The data is taken directly from the `reference` key within the provided node-data JSONs.
                These arguments are mutually exclusive.
                """
        ).add_mutually_exclusive_group()
    root_sequence = root_sequence_group.add_mutually_exclusive_group()
    root_sequence.add_argument('--include-root-sequence', action="store_true", help="Export as an additional JSON. The filename will follow the pattern of <OUTPUT>_root-sequence.json for a main auspice JSON of <OUTPUT>.json")
    root_sequence.add_argument('--include-root-sequence-inline', action="store_true", help="Export the root sequence within the dataset JSON. This should only be used for small genomes for file size reasons.")

    optional_settings = parser.add_argument_group(
        title="OTHER OPTIONAL SETTINGS"
    )
    optional_settings.add_argument(
        '--validation-mode',
        dest="validation_mode",
        type=ValidationMode,
        choices=[mode for mode in ValidationMode],
        default=ValidationMode.ERROR,
        help=validation_mode_help_message)
    optional_settings.add_argument(
        '--skip-validation',
        dest="validation_mode",
        action="store_const",
        const=ValidationMode.SKIP,
        help="Skip validation of input/output files, equivalent to --validation-mode=skip. Use at your own risk!")

    return parser


def set_display_defaults(data_json, config):
    # Note: these cannot be provided via command line args
    if config.get("display_defaults"):
        defaults = config["display_defaults"]
        if config.get("defaults"):
            deprecated("[config file] both 'defaults' (deprecated) and 'display_defaults' provided. Ignoring the former.")
    elif config.get("defaults"):
        deprecated("[config file] 'defaults' has been replaced with 'display_defaults'")
        defaults = config["defaults"]
    else:
        return

    v1_v2_keys = [ # each item: [0] v2 key name. [1] deprecated v1 key name
        ["geo_resolution", "geoResolution"],
        ["color_by", "colorBy"],
        ["distance_measure", "distanceMeasure"],
        ["map_triplicate", "mapTriplicate"]
    ]

    for [v2_key, v1_key] in [x for x in v1_v2_keys if x[1] in defaults]:
        deprecated("[config file] '{}' has been replaced with '{}'".format(v1_key, v2_key))
        defaults[v2_key] = defaults[v1_key]
        del defaults[v1_key]

    if defaults:
        data_json['meta']["display_defaults"] = defaults

def set_maintainers(data_json, config, cmd_line_maintainers):
    # Command-line args overwrite the config file
    # Command-line info could come in as multiple lists w/multiple values, ex:
    #       [['Name1 <url1>'], ['Name2 <url2>', 'Name3 <url3>'], ['Name4 <url4>']]
    # They may or may not all have URLs
    if cmd_line_maintainers:
        maintainers = []
        for arg_entry in cmd_line_maintainers:
            for maint in arg_entry:
                res = re.search('<(.*)>', maint)
                url = res.group(1) if res else ''
                name = maint.split("<")[0].strip()
                tmp_dict = {'name': name}
                if url:
                    tmp_dict['url'] = url
                maintainers.append(tmp_dict)
        data_json['meta']['maintainers'] = maintainers
    elif config.get("maintainer"): # v1-type specification
        data_json['meta']["maintainers"] = [{ "name": config["maintainer"][0], "url": config["maintainer"][1]}]
    elif config.get("maintainers"): # see schema for details
        data_json['meta']['maintainers'] = config['maintainers']
    else:
        warn("You didn't provide information on who is maintaining this analysis.")


def set_title(data_json, config, cmd_line_title):
    # title is not necessary. Cmd line args override any config settings
    if cmd_line_title:
        data_json['meta']['title'] = cmd_line_title
    elif config.get("title"):
        data_json['meta']['title'] = config.get("title")

def set_build_url(data_json, config, cmd_line_build_url):
    # build_url is not necessary. Cmd line args override any config settings
    if cmd_line_build_url:
        data_json['meta']['build_url'] = cmd_line_build_url
    elif config.get("build_url"):
        data_json['meta']['build_url'] = config.get("build_url")

def set_description(data_json, cmd_line_description_file):
    """
    Read Markdown file provided by *cmd_line_description_file* and set
    `meta.description` in *data_json* to the text provided.
    """
    try:
        with open_file(cmd_line_description_file) as description_file:
            markdown_text = description_file.read()
            data_json['meta']['description'] = markdown_text
    except FileNotFoundError:
        fatal("Provided desciption file {} does not exist".format(cmd_line_description_file))

def create_branch_mutations(branch_attrs, node_data):
    for node_name, node_info in node_data['nodes'].items():
        if node_name not in branch_attrs:
            continue # strain name not in the tree
        if "aa_muts" not in node_info and "muts" not in node_info:
            continue
        branch_attrs[node_name]['mutations'] = {}
        if "muts" in node_info and len(node_info["muts"]):
            branch_attrs[node_name]["mutations"]["nuc"] = node_info["muts"]
        if "aa_muts" in node_info:
            aa = {gene:data for gene, data in node_info["aa_muts"].items() if len(data)}
            branch_attrs[node_name]["mutations"].update(aa)

def create_branch_labels(branch_attrs, node_data, branch_data):
    ## start by creating the 'aa' branch label, summarising any amino acid mutations.
    ## (We have already set mutations on 'branch_attrs' if they exist, just not the label)
    ## This is done first so that if the user defines their own 'aa' labels they will
    ## overwrite the ones created here
    for branch_info in branch_attrs.values():
        genes = [gene for gene in branch_info.get('mutations', {}) if gene!='nuc']
        if len(genes):
            branch_info['labels']['aa'] = \
                '; '.join(f"{gene}: {', '.join(branch_info['mutations'][gene])}" for gene in genes)

    ## check for the special key 'clade_annotation' defined via node data.
    ## For historical reasons, this is interpreted as a branch label 'clade'
    for node_name, node_info in node_data.items():
        if node_name in branch_attrs and "clade_annotation" in node_info and is_valid(node_info["clade_annotation"]):
            branch_attrs[node_name]['labels']['clade'] = node_info["clade_annotation"]

    ## finally transfer any labels defined via <NODE DATA JSON> -> 'branches' -> labels
    for node_name, branch_info in branch_data.items():
        if node_name not in branch_attrs:
            continue
        for label_key, label_value in branch_info.get('labels', {}).items():
            if label_key.upper() == "NONE" or not is_valid(label_value):
                continue
            branch_attrs[node_name]["labels"][label_key] = label_value

def parse_node_data_and_metadata(T, node_data, metadata):
    node_data_names = set()
    metadata_names = set()

    # assign everything to node_attrs, exclusions considered later
    node_attrs = {clade.name: {} for clade in T.root.find_clades()}

    # first pass: metadata
    for metadata_id, node in metadata.items():
        if metadata_id in node_attrs: # i.e. this node name is in the tree
            for key, value in node.items():
                corrected_key = update_deprecated_names(key)
                node_attrs[metadata_id][corrected_key] = value
                metadata_names.add(corrected_key)

    # second pass: node data JSONs (overwrites keys of same name found in metadata)
    node_attrs_which_are_actually_branch_attrs = ["clade_annotation", "aa_muts", "muts"]
    for name, info in node_data['nodes'].items():
        if name in node_attrs: # i.e. this node name is in the tree
            for key, value in info.items():
                if key in node_attrs_which_are_actually_branch_attrs:
                    continue # these will be handled below
                corrected_key = update_deprecated_names(key)
                node_attrs[name][corrected_key] = value
                node_data_names.add(corrected_key)

    # third pass: create `branch_attrs`. The data comes from
    # (a) some keys within `node_data['nodes']` (for legacy reasons)
    # (b) the `node_data['branches']` dictionary, which currently only defines labels
    branch_attrs = {clade.name: defaultdict(dict) for clade in T.root.find_clades()}
    create_branch_mutations(branch_attrs, node_data)
    create_branch_labels(branch_attrs, node_data['nodes'], node_data.get('branches', {}))

    return (node_data, node_attrs, node_data_names, metadata_names, branch_attrs)

def get_config(args):
    if not args.auspice_config:
        return {}
    config = read_config(args.auspice_config)
    if args.validation_mode is not ValidationMode.SKIP:
        try:
            print("Validating config file {} against the JSON schema".format(args.auspice_config))
            validate_auspice_config_v2(args.auspice_config)
        except ValidateError:
            print("Validation of {} failed. Please check the formatting of this file & refer to the augur documentation for further help. ".format(args.auspice_config))
            validation_failure(args.validation_mode)
    # Print a warning about the inclusion of "vaccine_choices" which are _unused_ by `export v2`
    # (They are in the schema as this allows v1-compat configs to be used)
    if config.get("vaccine_choices"):
        warning("The config JSON can no longer specify the `vaccine_choices`, they must be specified through a node-data JSON. This info will be unused.")
        del config["vaccine_choices"]
    return config


def get_additional_metadata_columns(config, command_line_metadata_columns, metadata_names):
    # Command line args override what is set in the config file
    if command_line_metadata_columns:
        potential_metadata_columns = command_line_metadata_columns
    else:
        potential_metadata_columns = config.get("metadata_columns", [])

    additional_metadata_columns = []
    for col in potential_metadata_columns:
        # Match the column names corrected within parse_node_data_and_metadata
        corrected_col = update_deprecated_names(col)
        if corrected_col not in metadata_names:
            warn(f"Requested metadata column {col!r} does not exist and will not be exported")
            continue
        additional_metadata_columns.append(corrected_col)

    return additional_metadata_columns


def run(args):
    configure_warnings()
    data_json = {"version": "v2", "meta": {"updated": time.strftime('%Y-%m-%d')}}

    #load input files
    if args.node_data is not None:
      try:
          node_data_file = read_node_data(args.node_data, validation_mode=args.validation_mode) # node_data_files is an array of multiple files (or a single file)
      except FileNotFoundError:
          print(f"ERROR: node data file ({args.node_data}) does not exist")
          sys.exit(2)
    else:
        node_data_file = {'nodes': {}}

    if args.metadata is not None:
        try:
            metadata_df = read_metadata(
                args.metadata,
                delimiters=args.metadata_delimiters,
                id_columns=args.metadata_id_columns)

            # Add the index as a column.
            metadata_df[metadata_df.index.name] = metadata_df.index

            metadata_file = metadata_df.to_dict(orient="index")
        except FileNotFoundError:
            print(f"ERROR: meta data file ({args.metadata}) does not exist", file=sys.stderr)
            sys.exit(2)
        except InvalidDelimiter:
            raise AugurError(
                f"Could not determine the delimiter of {args.metadata!r}. "
                f"Valid delimiters are: {args.metadata_delimiters!r}. "
                "This can be changed with --metadata-delimiters."
            )
        except Exception as error:
            print(f"ERROR: {error}", file=sys.stderr)
            sys.exit(1)
    else:
        metadata_file = {}

    # parse input files
    T = Phylo.read(args.tree, 'newick')
    node_data, node_attrs, node_data_names, metadata_names, branch_attrs = \
            parse_node_data_and_metadata(T, node_data_file, metadata_file)
    config = get_config(args)
    additional_metadata_columns = get_additional_metadata_columns(config, args.metadata_columns, metadata_names)

    # set metadata data structures
    set_title(data_json, config, args.title)
    set_display_defaults(data_json, config)
    set_maintainers(data_json, config, args.maintainers)
    set_build_url(data_json, config, args.build_url)
    set_annotations(data_json, node_data)
    if args.description:
        set_description(data_json, args.description)

    try:
        set_colorings(
            data_json=data_json,
            config=get_config_colorings_as_dict(config),
            command_line_colorings=args.color_by_metadata,
            metadata_names=metadata_names,
            node_data_colorings=node_data_names,
            provided_colors=read_colors(args.colors),
            node_attrs=node_attrs,
            branch_attrs=branch_attrs
        )
    except FileNotFoundError as e:
        print(f"ERROR: required file could not be read: {e}")
        sys.exit(2)
    set_filters(data_json, config)

    # set tree structure
    data_json["tree"] = convert_tree_to_json_structure(T.root, node_attrs, node_div(T, node_attrs))
    set_node_attrs_on_tree(data_json, node_attrs, additional_metadata_columns)
    set_branch_attrs_on_tree(data_json, branch_attrs)

    set_geo_resolutions(data_json, config, args.geo_resolutions, read_lat_longs(args.lat_longs), node_attrs)
    set_panels(data_json, config, args.panels)
    set_data_provenance(data_json, config)

    # pass through any extensions block in the auspice config JSON without any changes / checking
    if config.get("extensions"):
        data_json["meta"]["extensions"] = config["extensions"]

    # Should output be minified?
    # User-specified arguments take precedence before determining behavior based
    # on the size of the tree.
    if args.minify_json or os.environ.get("AUGUR_MINIFY_JSON"):
        minify = True
    elif args.no_minify_json:
        minify = False
    else:
        if json_size(data_json) > MINIFY_THRESHOLD_MB * 10**6:
            minify = True
        else:
            minify = False

    # Write outputs - the (unified) dataset JSON intended for auspice & perhaps the ref root-sequence JSON
    indent = {"indent": None} if minify else {}
    if args.include_root_sequence or args.include_root_sequence_inline:
        # Note - argparse enforces that only one of these args will be true
        if 'reference' in node_data:
            if args.include_root_sequence_inline:
                data_json['root_sequence'] = node_data['reference']
            elif args.include_root_sequence:
                # Save the root sequence with the same stem from the main auspice
                # output filename.  For example, if the main auspice output is
                # "auspice/zika.json", the root sequence will be
                # "auspice/zika_root-sequence.json".
                output_path = Path(args.output)
                root_sequence_path = output_path.parent / Path(output_path.stem + "_root-sequence" + output_path.suffix)
                write_json(data=node_data['reference'], file=root_sequence_path, include_version=False, **indent)
        else:
            fatal("Root sequence output was requested, but the node data provided is missing a 'reference' key.")
    write_json(data=orderKeys(data_json), file=args.output, include_version=False, **indent)

    # validate outputs
    validate_data_json(args.output, args.validation_mode)

    if deprecationWarningsEmitted:
        print("\n------------------------")
        print("There were deprecation warnings displayed. They have been fixed but these will likely become breaking errors in a future version of augur.")
        print("------------------------")
    print("")
