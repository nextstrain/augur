"""
Functions that are from earlier versions of `augur validate` and haven't yet been
moved into `validate.py`. There's a lot of duplication here which can be abstracted
and refactored over time.
"""

import sys
from collections import defaultdict

def ensure_no_duplicate_names(root, ValidateError):
    """
    Check that all node names are identical, which is required for auspice (v2) JSONs.
    """
    names = set()
    def recurse(node):
        if node["name"] in names:
            raise ValidateError(f"Node {node['name']} appears multiple times in the tree.")
        names.add(node["name"])
        if "children" in node:
            [recurse(child) for child in node["children"]]
    recurse(root)


def collectTreeAttrsV2(root, warn):
    """
    Collect all keys specified on `node["node_attrs"]` throughout the tree
    and the values associated with them. Note that this will only look at
    attributes which are themselves objects with a `value` property.
    I.e. a node attribute `node["node_attrs"]["div"] -> numeric` will not
    be collected.
    Returns a tuple.
    return[0]: dict of `node_attr_property` -> x, where x is a dict with
    keys `count` -> INT, `values` -> SET, `onAllNodes` -> BOOL.
    return[1]: INT of number of terminal nodes in tree
    """
    seen = defaultdict(lambda: {"count": 0, "values": set(), "onAllNodes": False})
    num_nodes, num_terminal = (0, 0)
    def recurse(node):
        nonlocal num_nodes, num_terminal
        num_nodes += 1
        for prop, info in node.get("node_attrs", {}).items():
            if not isinstance(info, dict) or "value" not in info:
                continue
            seen[prop]["count"] += 1
            seen[prop]["values"].add(info["value"])
        if "children" in node:
            [recurse(child) for child in node["children"]]
        else:
            num_terminal += 1
    recurse(root)

    for data in seen.values():
        if data["count"] == num_nodes:
            data["onAllNodes"] = True

    return(seen, num_terminal)


def collectMutationGenes(root):
    """
    Returns a set of all genes specified in the tree in the "aa_muts" objects
    """
    genes = set()
    def recurse(node):
        mutations = node.get("branch_attrs", {}).get("mutations", False)
        if mutations:
            genes.update(mutations.keys())
        if "children" in node:
            [recurse(child) for child in node["children"]]
    recurse(root)
    genes -= {"nuc"}
    return genes

def collectBranchLabels(root):
    labels = set()
    def recurse(node):
        labels.update(node.get("branch_attrs", {}).get("labels", {}).keys())
        if "children" in node:
            [recurse(child) for child in node["children"]]
    recurse(root)
    return labels

def verifyMainJSONIsInternallyConsistent(data, ValidateError):
    """
    Check possible sources of conflict within the main (unified) JSON
    This function is only used for schema v2.0
    In an ideal world, the JSON schema would be able to validate this,
    however this function performs tests which are not possible to
    define in the schema (as of JSON schema v6).
    """
    warnings = False
    def warn(msg):
        nonlocal warnings
        warnings = True
        print("\tWARNING: ", msg, file=sys.stderr)

    print("Validating that the JSON is internally consistent...")

    ensure_no_duplicate_names(data["tree"], ValidateError)

    if "entropy" in data["meta"]["panels"] and "genome_annotations" not in data["meta"]:
        warn("The entropy panel has been specified but annotations don't exist.")

    tree_traits, _ = collectTreeAttrsV2(data["tree"], warn)

    if "geo_resolutions" in data["meta"]:
        for geo_res in data["meta"]["geo_resolutions"]:
            geo_name = geo_res["key"]
            deme_to_lat_longs = geo_res["demes"]
            if geo_name not in tree_traits:
                warn("The geographic resolution \"{}\" does not appear on any tree nodes.".format(geo_name))
                continue
            # pass 1: check the demes in "geo_resolutions" are found on the tree
            for geoValue in deme_to_lat_longs.keys():
                if geoValue not in tree_traits[geo_name]["values"]:
                    warn("\"{}\", a value of the geographic resolution \"{}\", does not appear on any tree nodes.".format(geoValue, geo_name))
            # pass 1: check the demes across the tree are represented in "geo_resolutions"
            for geoValue in tree_traits[geo_name]["values"]:
                if geoValue not in deme_to_lat_longs:
                    warn("\"{}\", a value of the geographic resolution \"{}\", appears in the tree but not in the metadata."
                        "\n\t\tThis will cause transmissions & demes involving this location not to be displayed in Auspice".format(geoValue, geo_name))
    else:
        if "map" in data["meta"]["panels"]:
            warn("Map panel was requested but no geographic_info was provided")

    if "colorings" in data["meta"]:
        for coloring in data["meta"]["colorings"]:
            colorBy = coloring["key"]
            if colorBy == "gt":
                continue
            if colorBy not in tree_traits:
                warn("The coloring \"{}\" does not appear as an attr on any tree nodes.".format(colorBy))
            if "scale" in coloring:
                scale = coloring["scale"]
                if isinstance(scale, list):
                    for value, hex in scale:
                        if value not in tree_traits[colorBy]["values"]:
                            warn("Color option \"{}\" specifies a hex code for \"{}\" but this isn't ever seen on the tree nodes.".format(colorBy, value))
                elif isinstance(scale, str):
                    raise ValidateError("String colour scales are not yet implemented")
                else:
                    raise ValidateError("Invalid color scale (for trait \"{}\")".format(colorBy))
            if "domain" in coloring:
                domain = coloring["domain"]
                if coloring["type"] in ["ordinal", "categorical"]:
                    inMetaNotInTree = [val for val in domain if val not in tree_traits[colorBy]["values"]]
                    if len(inMetaNotInTree):
                        warn("Domain for {} defined the following values which are not present on the tree: {}".format(colorBy, inMetaNotInTree.join(", ")))
                    inTreeNotInMeta = [val for val in tree_traits[colorBy]["values"] if val not in domain]
                    if len(inTreeNotInMeta):
                        warn("Tree defined values for {} which were not in the domain: {}".format(colorBy, inTreeNotInMeta.join(", ")))
                elif coloring["type"] == "boolean":
                    raise ValidateError("Cannot povide a domain for a boolean coloring ({})".format(colorBy))
    else:
        warn("No colourings were provided")

    if "filters" in data["meta"]:
        for filter in data["meta"]["filters"]:
            if filter not in tree_traits:
                warn("The filter \"{}\" does not appear as a property on any tree nodes.".format(filter))

    genes_with_mutations = collectMutationGenes(data['tree'])
    if len(genes_with_mutations):
        if "genome_annotations" not in data["meta"]:
            warn("The tree defined mutations on genes {}, but annotations aren't defined in the meta JSON.".format(", ".join(genes_with_mutations)))
        else:
            for gene in genes_with_mutations:
                if gene not in data["meta"]["genome_annotations"]:
                    warn("The tree defined mutations on gene {} which doesn't appear in the metadata annotations object.".format(gene))

    default_branch_label = data.get("meta").get("display_defaults", {}).get("branch_label")
    if default_branch_label and default_branch_label.lower() != "none":
        labels = collectBranchLabels(data['tree'])
        if not default_branch_label in labels:
            warn("Default label to display \"{}\" isn't found anywhere on the tree!".format(default_branch_label))

    return not warnings


def collectTreeAttrsV1(root):
    """
    Collect all keys specified on node->attr (or node->traits) throughout the tree
    If the values of these keys are strings, then also collect the values
    """
    seen = defaultdict(lambda: {"count": 0, "values": set(), "onAllNodes": False})
    num_nodes, num_terminal = (0, 0)
    def recurse(node):
        nonlocal num_nodes, num_terminal
        num_nodes += 1
        traits = node["attr"].keys()
        for property in traits:
            seen[property]["count"] += 1
            value = node["attr"][property]
            if isinstance(value, str):
                seen[property]["values"].add(value)

        if "children" in node:
            [recurse(child) for child in node["children"]]
        else:
            num_terminal += 1
    recurse(root)

    for data in seen.values():
        if data["count"] == num_nodes:
            data["onAllNodes"] = True

    return(seen, num_terminal)


def collectAAMutationGenesV1(root):
    """
    Returns a set of all genes specified in the tree in the "aa_muts" objects
    """
    genes = set()
    def recurse(node):
        if "aa_muts" in node and len(node["aa_muts"]):
            genes.update(node["aa_muts"].keys())
        if "children" in node:
            [recurse(child) for child in node["children"]]
    recurse(root)
    return genes


def verifyMetaAndOrTreeJSONsAreInternallyConsistent(meta_json, tree_json, ValidateError):
    """
    Check all possible sources of conflict internally & between the metadata & tree JSONs
    This is only that which cannot be checked by the schemas
    """
    warnings = False
    def warn(msg):
        nonlocal warnings
        warnings = True
        print("\tWARNING: ", msg, file=sys.stderr)


    print("Validating that meta + tree JSONs are internally consistent...")
    mj = meta_json

    if "panels" in mj and "entropy" in mj["panels"] and "annotations" not in mj:
        warn("\tERROR: The entropy panel has been specified but annotations don't exist.")

    tj = tree_json
    treeAttrs, num_terminal_nodes = collectTreeAttrsV1(tj)

    if "geo" in mj:
        for geoName in mj["geo"]:
            if geoName not in treeAttrs:
                warn("The geographic resolution \"{}\" does not appear as an attr on any tree nodes.".format(geoName))
            else:
                for geoValue in mj["geo"][geoName]:
                    if geoValue not in treeAttrs[geoName]["values"]:
                        warn("\"{}\", a value of the geographic resolution \"{}\", does not appear as a value of attr->{} on any tree nodes.".format(geoValue, geoName, geoName))
                for geoValue in treeAttrs[geoName]["values"]:
                    if geoValue not in mj["geo"][geoName]:
                        warn("\"{}\", a value of the geographic resolution \"{}\", appears in the tree but not in the metadata.".format(geoValue, geoName))
                        warn("\tThis will cause transmissions & demes involving this location not to be displayed in Auspice")


    if "color_options" in mj:
        for colorBy in [x for x in mj["color_options"] if x != "gt"]:
            if colorBy not in treeAttrs:
                warn("The color_option \"{}\" does not appear as an attr on any tree nodes.".format(colorBy))
            elif "color_map" in mj["color_options"][colorBy]:
                # if there's a color_map, then check that there are no values which aren't seen on the tree
                for (value, hex) in mj["color_options"][colorBy]["color_map"]:
                    if value not in treeAttrs[colorBy]["values"]:
                        warn("Color option \"{}\" specifies a hex code for \"{}\" but this isn't ever seen on the tree nodes.".format(colorBy, value))
                # inversely, check for values on the tree not defined in the color_map
                for value in treeAttrs[colorBy]["values"]:
                    color_map_values = [x[0] for x in mj["color_options"][colorBy]["color_map"]]
                    if value not in color_map_values:
                        warn("Color option \"{}\", which contains a color_map, is missing \"{}\"".format(colorBy, value))



    if "filters" in mj:
        for filter in mj["filters"]:
            if filter not in treeAttrs:
                warn("The filter \"{}\" does not appear on any tree nodes.".format(filter))


    if "virus_count" in mj and mj["virus_count"] != num_terminal_nodes:
        raise ValidateError("Meta JSON virus_count ({}) differs from the number of nodes in the tree ({})".format(mj["virus_count"], num_terminal_nodes))

    genes_with_aa_muts = collectAAMutationGenesV1(tj)
    if len(genes_with_aa_muts):
        if "annotations" not in mj:
            warn("The tree defined AA mutations on genes {}, but annotations aren't defined in the meta JSON.".format(", ".join(genes_with_aa_muts)))
        else:
            for gene in genes_with_aa_muts:
                if gene not in mj["annotations"]:
                    warn("The tree defined AA mutations on gene {} which doesn't appear in the metadata annotations object.".format(gene))

    return not warnings
