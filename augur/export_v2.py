"""
Export JSON files suitable for visualization with auspice.
"""
from pathlib import Path
import os, sys
import time
from collections import defaultdict, deque
import warnings
import numbers
import re
from Bio import Phylo
from .utils import read_metadata, read_node_data, write_json, read_config, read_lat_longs, read_colors
from .validate import export_v2 as validate_v2, auspice_config_v2 as validate_auspice_config_v2, ValidateError

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

def convert_tree_to_json_structure(node, metadata, div=0):
    """
    converts the Biopython tree structure to a dictionary that can
    be written to file as a json. This is called recursively.
    Creates the name property & divergence on each node

    input
        node -- node for which top level dict is produced.
        div  -- cumulative divergence (root = 0). False → divergence won't be exported.

    returns
        tree in JSON structure
        list of strains
    """

    # Does the tree have divergence? (BEAST trees may not)
    # only calculate this for the root node!
    if div == 0 and 'mutation_length' not in metadata[node.name] and 'branch_length' not in metadata[node.name]:
        div = False

    node_struct = {'name': node.name, 'node_attrs': {}, 'branch_attrs': {}}
    if div is not False: # div=0 is ok
        node_struct["node_attrs"]["div"] = div

    if node.clades:
        node_struct["children"] = []
        for child in node.clades:
            if div is False:
                cdiv=False
            else:
                if 'mutation_length' in metadata[child.name]:
                    cdiv = div + metadata[child.name]['mutation_length']
                elif 'branch_length' in metadata[child.name]:
                    cdiv = div + metadata[child.name]['branch_length']
                else:
                    print("ERROR: Cannot find branch length information for %s"%(child.name))

            node_struct["children"].append(convert_tree_to_json_structure(child, metadata, div=cdiv))

    return node_struct

def are_mutations_defined(node_attrs):
    for node, data in node_attrs.items():
        if data.get("aa_muts") or data.get("muts"):
            return True
    return False


def are_clades_defined(node_attrs):
    for node, data in node_attrs.items():
        if data.get("clade_membership") or data.get("clade_annotation"):
            return True
    return False


def are_dates_defined(node_attrs):
    for node, data in node_attrs.items():
        if data.get("num_date"):
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
        if data.get(key):
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


def set_colorings(data_json, config, command_line_colorings, metadata_names, node_data_colorings, provided_colors, node_attrs):

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
            allowedTypes = ["continuous", "ordinal", "categorical", "boolean"]
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
        if scale_type is "continuous":
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

    def _create_coloring(key):
        # handle deprecations
        if key == "authors":
            deprecated("[colorings] The 'authors' key is now called 'author'")
            key = "author"
        return {"key": key}

    def _is_valid(coloring):
        key = coloring["key"]
        trait_values = get_values_across_nodes(node_attrs, key) # e.g. list of countries, regions etc
        if key == "gt" and not are_mutations_defined(node_attrs):
            warn("[colorings] You asked for mutations (\"gt\"), but none are defined on the tree. They cannot be used as a coloring.")
            return False
        if key != "gt" and not trait_values:
            warn("You asked for a color-by for trait '{}', but it has no values on the tree. It has been ignored.".format(key))
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
                colorings.append(_create_coloring(x))
            # then add in command line colorings
            for x in command_line_colorings:
                colorings.append(_create_coloring(x))
        else:
            # if we have a config file, start with these (extra info, such as title&type, is added in later)
            if config:
                for x in config.keys():
                    colorings.append(_create_coloring(x))
            # then add in any auto-colorings already validated to include
            for x in auto_colorings:
                colorings.append(_create_coloring(x))

        explicitly_defined_colorings = [x["key"] for x in colorings]
        # add in genotype as a special case if (a) not already set and (b) the data supports it
        if "gt" not in explicitly_defined_colorings and are_mutations_defined(node_attrs):
            colorings.insert(0,{'key':'gt'})
        if "num_date" not in explicitly_defined_colorings and are_dates_defined(node_attrs):
            colorings.insert(0,{'key':'num_date'})
        if "clade_membership" not in explicitly_defined_colorings and are_clades_defined(node_attrs):
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

def validate_data_json(filename):
    print("Validating produced JSON")
    try:
        validate_v2(main_json=filename)
    except ValidateError as e:
        print(e)
        print("\n------------------------")
        print("Validation of {} failed. Please check this in a local instance of `auspice`, as it is not expected to display correctly. ".format(filename))
        print("------------------------")
        sys.exit(2)


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


def set_node_attrs_on_tree(data_json, node_attrs):
    '''
    Assign desired colorings, metadata etc to the tree structure

    Parameters
    ----------
    data_json : dict
    node_attrs: dict
        keys: strain names. values: dict with keys -> all available metadata (even "excluded" keys), values -> data (string / numeric / bool)
    '''

    author_data = create_author_data(node_attrs)

    def _transfer_mutations(node, raw_data):
        if "aa_muts" in raw_data or "muts" in raw_data:
            node["branch_attrs"]["mutations"] = {}
            if "muts" in raw_data and len(raw_data["muts"]):
                node["branch_attrs"]["mutations"]["nuc"] = raw_data["muts"]
            if "aa_muts" in raw_data:
                aa = {gene:data for gene, data in raw_data["aa_muts"].items() if len(data)}
                node["branch_attrs"]["mutations"].update(aa)
                #convert mutations into a label
                if aa:
                    aa_lab = '; '.join("{!s}: {!s}".format(key,', '.join(val)) for (key,val) in aa.items())
                    if 'labels' in node["branch_attrs"]:
                        node["branch_attrs"]["labels"]["aa"] = aa_lab
                    else:
                        node["branch_attrs"]["labels"] = { "aa": aa_lab }

    def _transfer_vaccine_info(node, raw_data):
        if raw_data.get("vaccine"):
            node["node_attrs"]['vaccine'] = raw_data['vaccine']

    def _transfer_labels(node, raw_data):
        if "clade_annotation" in raw_data and is_valid(raw_data["clade_annotation"]):
            if 'labels' in node["branch_attrs"]:
                node["branch_attrs"]["labels"]['clade'] = raw_data["clade_annotation"]
            else:
                node["branch_attrs"]["labels"] = { "clade": raw_data["clade_annotation"] }

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
            node["node_attrs"]["num_date"] = {"value": raw_data["num_date"]}
            if is_valid(raw_data.get("num_date_confidence", None)):
                node["node_attrs"]["num_date"]["confidence"] = raw_data["num_date_confidence"]

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
            if is_valid(raw_data.get(key, None)):
                node["node_attrs"][key] = {"value": raw_data[key]}
                if is_valid(raw_data.get(key+"_confidence", None)):
                    node["node_attrs"][key]["confidence"] = raw_data[key+"_confidence"]
                if is_valid(raw_data.get(key+"_entropy", None)):
                    node["node_attrs"][key]["entropy"] = raw_data[key+"_entropy"]

    def _transfer_author_data(node):
        if node["name"] in author_data:
            node["node_attrs"]["author"] = author_data[node["name"]]

    def _recursively_set_data(node):
        # get all the available information for this particular node
        raw_data = node_attrs[node["name"]]
        # transfer "special cases"
        _transfer_mutations(node, raw_data)
        _transfer_vaccine_info(node, raw_data)
        _transfer_labels(node, raw_data)
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
        "clade_annotation", # Clade annotation is label, not colorby!
        "clade_membership", # will be auto-detected if it is available
        "authors",          # authors are set as a node property, not a trait property
        "author",           # see above
        "vaccine",          # vaccine info is stored as a "special" node prop
        'branch_length',
        'num_date',
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

def get_root_sequence(root_node, ref=None, translations=None):
    '''
    create a json structure that contains the sequence of the root, both as
    nucleotide and as translations. This allows look-up of the sequence for
    all states, including those that are not variable.

    Parameters
    ----------
    root_node : dict
    	data associated with the node
    ref : str, optional
        filename of the root sequence
    translations : str, optional
        file name of translations

    Returns
    -------
    dict
        dict of nucleotide sequence and translations
    '''
    root_sequence = {}
    if ref and translations:
        from Bio import SeqIO
        refseq = SeqIO.read(ref, 'fasta')
        root_sequence['nuc']=str(refseq.seq)
        for gene in SeqIO.parse(translations, 'fasta'):
            root_sequence[gene.id] = str(gene.seq)
    else:
        root_sequence["nuc"] = root_node["sequence"]
        root_sequence.update(root_node["aa_sequences"])

    return root_sequence


def register_arguments_v2(subparsers):
    v2 = subparsers.add_parser("v2", help="Export version 2 JSON schema")

    required = v2.add_argument_group(
        title="REQUIRED"
    )
    required.add_argument('--tree','-t', metavar="newick", required=True, help="Phylogenetic tree, usually output from `augur refine`")
    required.add_argument('--node-data', metavar="JSON", required=True, nargs='+', help="JSON files containing metadata for nodes in the tree")
    required.add_argument('--output', metavar="JSON", required=True, help="Ouput file (typically for visualisation in auspice)")

    config = v2.add_argument_group(
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
    config.add_argument('--geo-resolutions', metavar="trait", nargs='+', help="Geographic traits to be displayed on map")
    config.add_argument('--color-by-metadata', metavar="trait", nargs='+', help="Metadata columns to include as coloring options")
    config.add_argument('--panels', metavar="panels", nargs='+', choices=['tree', 'map', 'entropy', 'frequencies'], help="Restrict panel display in auspice. Options are %(choices)s. Ignore this option to display all available panels.")

    optional_inputs = v2.add_argument_group(
        title="OPTIONAL INPUT FILES"
    )
    optional_inputs.add_argument('--metadata', metavar="FILE", help="Additional metadata for strains in the tree, as CSV or TSV")
    optional_inputs.add_argument('--colors', metavar="FILE", help="Custom color definitions, one per line in the format `TRAIT_TYPE\\tTRAIT_VALUE\\tHEX_CODE`")
    optional_inputs.add_argument('--lat-longs', metavar="TSV", help="Latitudes and longitudes for geography traits (overrides built in mappings)")

    optional_settings = v2.add_argument_group(
        title="OPTIONAL SETTINGS"
    )
    optional_settings.add_argument('--minify-json', action="store_true", help="export JSONs without indentation or line returns")
    optional_settings.add_argument('--include-root-sequence', action="store_true", help="Export an additional JSON containing the root sequence (reference sequence for vcf) used to identify mutations. The filename will follow the pattern of <OUTPUT>_root-sequence.json for a main auspice JSON of <OUTPUT>.json")

    return v2


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
        with open(cmd_line_description_file, encoding='utf-8') as description_file:
            markdown_text = description_file.read()
            data_json['meta']['description'] = markdown_text
    except FileNotFoundError:
        fatal("Provided desciption file {} does not exist".format(cmd_line_description_file))

def parse_node_data_and_metadata(T, node_data, metadata):
    node_data_names = set()
    metadata_names = set()

    # assign everything to node_attrs, exclusions considered later
    node_attrs = {clade.name: {} for clade in T.root.find_clades()}

    # first pass: metadata
    for node in metadata.values():
        if node["strain"] in node_attrs: # i.e. this node name is in the tree
            for key, value in node.items():
                corrected_key = update_deprecated_names(key)
                node_attrs[node["strain"]][corrected_key] = value
                metadata_names.add(corrected_key)

    # second pass: node data JSONs (overwrites keys of same name found in metadata)
    for name, info in node_data['nodes'].items():
        if name in node_attrs: # i.e. this node name is in the tree
            for key, value in info.items():
                corrected_key = update_deprecated_names(key)
                node_attrs[name][corrected_key] = value
                node_data_names.add(corrected_key)

    return (node_data, node_attrs, node_data_names, metadata_names)

def get_config(args):
    if not args.auspice_config:
        return {}
    config = read_config(args.auspice_config)
    try:
        print("Validating config file {} against the JSON schema".format(args.auspice_config))
        validate_auspice_config_v2(args.auspice_config)
    except ValidateError:
        print("Validation of {} failed. Please check the formatting of this file & refer to the augur documentation for further help. ".format(args.auspice_config))
        sys.exit(2)
    # Print a warning about the inclusion of "vaccine_choices" which are _unused_ by `export v2`
    # (They are in the schema as this allows v1-compat configs to be used)
    if config.get("vaccine_choices"):
        warning("The config JSON can no longer specify the `vaccine_choices`, they must be specified through a node-data JSON. This info will be unused.")
        del config["vaccine_choices"]
    return config

def run_v2(args):
    configure_warnings()
    data_json = {"version": "v2", "meta": {"updated": time.strftime('%Y-%m-%d')}}

    #load input files
    try:
        node_data_file = read_node_data(args.node_data) # node_data_files is an array of multiple files (or a single file)
    except FileNotFoundError:
        print(f"ERROR: node data file ({args.node_data}) does not exist")
        sys.exit(2)

    if args.metadata is not None:
        try:
            metadata_file, _ = read_metadata(args.metadata)
        except FileNotFoundError:
            print(f"ERROR: meta data file ({args.metadata}) does not exist")
            sys.exit(2)
    else:
        metadata_file = {}

    # parse input files
    T = Phylo.read(args.tree, 'newick')
    node_data, node_attrs, node_data_names, metadata_names = parse_node_data_and_metadata(T, node_data_file, metadata_file)
    config = get_config(args)

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
            node_attrs=node_attrs
        )
    except FileNotFoundError as e:
        print(f"ERROR: required file could not be read: {e}")
        sys.exit(2)
    set_filters(data_json, config)

    # set tree structure
    data_json["tree"] = convert_tree_to_json_structure(T.root, node_attrs)
    set_node_attrs_on_tree(data_json, node_attrs)
    set_geo_resolutions(data_json, config, args.geo_resolutions, read_lat_longs(args.lat_longs), node_attrs)
    set_panels(data_json, config, args.panels)
    set_data_provenance(data_json, config)

    # Write outputs - the (unified) dataset JSON intended for auspice & perhaps the ref root-sequence JSON
    indent = {"indent": None} if args.minify_json else {}
    write_json(data=data_json, file_name=args.output, include_version=False, **indent)

    if args.include_root_sequence:
        if 'reference' in node_data:
            # Save the root sequence with the same stem from the main auspice
            # output filename.  For example, if the main auspice output is
            # "auspice/zika.json", the root sequence will be
            # "auspice/zika_root-sequence.json".
            output_path = Path(args.output)
            root_sequence_path = output_path.parent / Path(output_path.stem + "_root-sequence" + output_path.suffix)
            write_json(data=node_data['reference'], file_name=root_sequence_path, include_version=False, **indent)
        else:
            fatal("Root sequence output was requested, but the node data provided is missing a 'reference' key.")

    # validate outputs
    validate_data_json(args.output)

    if deprecationWarningsEmitted:
        print("\n------------------------")
        print("There were deprecation warnings displayed. They have been fixed but these will likely become breaking errors in a future version of augur.")
        print("------------------------")
    print("")
