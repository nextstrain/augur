"""
Validate a set of JSON files intended for visualization in auspice.
"""

import jsonschema
import os, sys
import json
from argparse import SUPPRESS
from collections import defaultdict
from pkg_resources import resource_string, resource_filename
from copy import deepcopy

class ValidateError(Exception):
    pass


def loadJSONsToValidate(paths):
    """
    paths: array of paths to JSONs.
    JSON type (and therefore schema type) is inferred from the pathname, paralleling Auspice
    tree & meta JSONs
    """
    ret = {}
    for path in paths:
        with open(path) as f:
            try:
                jsonToValidate = json.load(f)
            except json.JSONDecodeError:
                raise ValidateError("Supplied JSON to validate ({}) is not a valid JSON".format(path))
        if path.endswith("_tree.json"):
            schema_type = "tree"
        elif path.endswith("_meta.json"):
            schema_type = "meta"
        elif path.endswith("frequencies.json") or path.endswith("entropy.json") or path.endswith("sequences.json"):
            print("JSON to validate ({}) has no schema (yet). Continuing.".format(path))
            continue

        if schema_type in ret:
            raise ValidateError('Cannot deal with multiple JSONs of the same type ({})'.format(schema_type))

        ret[schema_type] = {
            "path": path,
            "json": jsonToValidate
        }
    return ret


def checkSchemaIsValid(schema, name):
    # see http://python-jsonschema.readthedocs.io/en/latest/errors/
    try:
        jsonschema.Draft6Validator.check_schema(schema)
    except jsonschema.exceptions.SchemaError as err:
        raise ValidateError("Schema {} is not a valid JSON file. Error: {}".format(path, err))


def loadSchemas(types):
    """
    For types (such as "tree", "meta"), load and internally-verify the schema.
    """
    ret = {}

    for schema_type in types:
        if schema_type in ret:
            continue
        if schema_type == "meta":
            path = "data/schema_meta.json"
        elif schema_type == "tree":
            path = "data/schema_tree.json"
        else:
            raise ValidateError("ERROR: No specified (v1) schema for ", schema_type)

        try:
            schema = json.loads(resource_string(__package__, path))
        except json.JSONDecodeError as err:
            raise ValidateError("Schema {} is not a valid JSON file. Error: {}".format(path, err))
        checkSchemaIsValid(schema, path)
        ret[schema_type] = jsonschema.Draft6Validator(schema)

    return ret


def collectTreeAttrs(root):
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


def collectAAMutationGenes(root):
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


def verifyMetaAndOrTreeJSONsAreInternallyConsistent(meta, tree):
    """
    Check all possible sources of conflict internally & between the metadata & tree JSONs
    This is only that which cannot be checked by the schemas
    """
    return_status = 0

    def error(msg):
        nonlocal return_status
        return_status = 1
        print("\tERROR: ", msg, file=sys.stderr)

    def warn(msg):
        print("\tWARNING: ", msg, file=sys.stderr)


    print("Validating that {} and {} are internally consistent... ".format(meta["path"], tree["path"]))
    mj = meta["json"]

    if "panels" in mj and "entropy" in mj["panels"] and "annotations" not in mj:
        error("\tERROR: The entropy panel has been specified but annotations don't exist.")

    if not tree:
        return return_status

    tj = tree["json"]
    treeAttrs, num_terminal_nodes = collectTreeAttrs(tj)

    if "geo" in mj:
        for geoName in mj["geo"]:
            if geoName not in treeAttrs:
                error("The geographic resolution \"{}\" does not appear as an attr on any tree nodes.".format(geoName))
            else:
                for geoValue in mj["geo"][geoName]:
                    if geoValue not in treeAttrs[geoName]["values"]:
                        warn("\"{}\", a value of the geographic resolution \"{}\", does not appear as a value of attr->{} on any tree nodes.".format(geoValue, geoName, geoName))
                for geoValue in treeAttrs[geoName]["values"]:
                    if geoValue not in mj["geo"][geoName]:
                        error("\"{}\", a value of the geographic resolution \"{}\", appears in the tree but not in the metadata.".format(geoValue, geoName))
                        error("\tThis will cause transmissions & demes involving this location not to be displayed in Auspice")


    if "color_options" in mj:
        for colorBy in [x for x in mj["color_options"] if x != "gt"]:
            if colorBy not in treeAttrs:
                error("The color_option \"{}\" does not appear as an attr on any tree nodes.".format(colorBy))
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
                error("The filter \"{}\" does not appear on any tree nodes.".format(filter))

    if "author_info" in mj:
        if "authors" not in treeAttrs:
            warn("\"author_info\" exists in the metadata JSON but \"authors\" are never defined on a tree node attr.")
        else:
            for author in mj["author_info"]:
                if author not in treeAttrs["authors"]["values"]:
                    warn("Author \"{}\" defined in \"author_info\" but does not appear on any tree node.".format(author))
            for author in treeAttrs["authors"]["values"]:
                if author not in mj["author_info"]:
                    error("Author \"{}\" defined on a tree node but not in \"author_info\".".format(author))

    if "virus_count" in mj and mj["virus_count"] != num_terminal_nodes:
        error("Meta JSON virus_count ({}) differs from the number of nodes in the tree ({})".format(mj["virus_count"], num_terminal_nodes))

    genes_with_aa_muts = collectAAMutationGenes(tj)
    if len(genes_with_aa_muts):
        if "annotations" not in mj:
            error("The tree defined AA mutations on genes {}, but annotations aren't defined in the meta JSON.".format(", ".join(genes_with_aa_muts)))
        else:
            for gene in genes_with_aa_muts:
                if gene not in mj["annotations"]:
                    error("The tree defined AA mutations on gene {} which doesn't appear in the metadata annotations object.".format(gene))

    return return_status


def register_arguments_v1(subparsers):
    v1 = subparsers.add_parser("v1", help="Validate version 1 JSONs")
    v1.add_argument('--json', required=True, nargs='+', help="JSONs to validate")
    v1.add_argument('--v1', default=True, help=SUPPRESS)


def run_v1(args):
    '''
    Validate auspice-compatable JSONs against a schema
    '''
    return_status = 0
    try:
        datasets = loadJSONsToValidate(args.json)
        schemas = loadSchemas([x for x in datasets])
    except ValidateError as e:
        print(e)
        return 1


    for schema_type, data in datasets.items():
        print("Validating {} using the {} schema (version 1)... ".format(data["path"], schema_type), end='')
        validator = schemas[schema_type]

        try:
            validator.validate(data["json"]) # https://python-jsonschema.readthedocs.io/en/latest/validate/
        except jsonschema.exceptions.ValidationError:
            print("FAILED")
            return_status = 1
            for error in sorted(validator.iter_errors(data["json"]), key=str):
                trace = list(error.schema_path) # this may be really long for nested tree structures
                if len(trace) > 6:
                    trace = ["..."] + trace[-5:]
                trace = [str(x) for x in trace]
                print("\tERROR: {}. Trace: {}".format(error.message, " - ".join(trace)), file=sys.stderr)
        else:
            print("SUCCESS")

    if "meta" in datasets and "tree" in datasets:
        return_status = verifyMetaAndOrTreeJSONsAreInternallyConsistent(datasets["meta"], datasets["tree"]) | return_status

    if not return_status:
        print("SUCCESS")
    return return_status
