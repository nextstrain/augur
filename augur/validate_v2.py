"""
Validate a set of JSON files intended for visualization in auspice.
"""

import jsonschema
import os, sys
import json
from collections import defaultdict
from pkg_resources import resource_string, resource_filename
from copy import deepcopy


class ValidateError(Exception):
    pass

def loadJSONsToValidate(path):
    """
    path: path to JSON.
    JSON type (and therefore schema type) is inferred from the pathname
    """
    ret = {}

    with open(path) as f:
        try:
            jsonToValidate = json.load(f)
        except json.JSONDecodeError:
            raise ValidateError("Supplied JSON to validate ({}) is not a valid JSON".format(path))

    if path.endswith("frequencies.json") or path.endswith("entropy.json") or path.endswith("sequences.json"):
        print("JSON to validate ({}) has no schema (yet).".format(path))
    else:
        schema_type = "main"

    if schema_type in ret:
        raise ValidateError('Cannot deal with multiple JSONs')

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
        if schema_type == "main":
            path = "data/schema.json"
        else:
            raise ValidateError("ERROR: No specified schema for ", schema_type)

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
        traits = node["traits"].keys()
        for property in traits:
            # Process author info from node not traits
            if property == "authors":
                continue
            seen[property]["count"] += 1
            value = node["traits"][property]["value"]
            if isinstance(value, str):
                seen[property]["values"].add(value)

        if "authors" in node:
            seen["authors"]["count"] += 1
            seen["authors"]["values"].add(node["authors"].lower())
        if "num_date" in node:
            seen["num_date"]["count"] += 1
            seen["num_date"]["values"].add(node["num_date"]["value"])

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
        if "mutations" in node:
            genes.update(node["mutations"].keys())
        if "children" in node:
            [recurse(child) for child in node["children"]]
    recurse(root)
    genes -= {"nuc"}
    return genes


def verifyMainJSONIsInternallyConsistent(main):
    """
    Check all possible sources of conflict within the main (unified) JSON
    This function is only used for schema v2.0
    """
    return_status = 0

    def error(msg):
        nonlocal return_status
        return_status = 1
        print("\tERROR: ", msg, file=sys.stderr)

    def warn(msg):
        print("\tWARNING: ", msg, file=sys.stderr)


    print("Validating that {} is internally consistent... ".format(main["path"]))
    data = main["json"]

    if "entropy" in data["panels"] and "genome_annotations" not in data:
        error("The entropy panel has been specified but annotations don't exist.")

    treeTraits, _ = collectTreeAttrs(data["tree"])

    if "geographic_info" in data:
        for geoName in data["geographic_info"]:
            if geoName not in treeTraits:
                error("The geographic resolution \"{}\" does not appear as an attr on any tree nodes.".format(geoName))
            else:
                for geoValue in data["geographic_info"][geoName]:
                    if geoValue not in treeTraits[geoName]["values"]:
                        warn("\"{}\", a value of the geographic resolution \"{}\", does not appear as a value of attr->{} on any tree nodes.".format(geoValue, geoName, geoName))
                for geoValue in treeTraits[geoName]["values"]:
                    if geoValue not in data["geographic_info"][geoName]:
                        error("\"{}\", a value of the geographic resolution \"{}\", appears in the tree but not in the metadata.".format(geoValue, geoName))
                        error("\tThis will cause transmissions & demes involving this location not to be displayed in Auspice")
    else:
        if "map" in data["panels"]:
            error("Map panel was requested but no geographic_info was provided")


    if "colorings" in data:
        for colorBy in [x for x in data["colorings"] if x != "gt"]:
            if colorBy not in treeTraits:
                error("The coloring \"{}\" does not appear as an attr on any tree nodes.".format(colorBy))
            if "scale" in data["colorings"][colorBy]:
                scale = data["colorings"][colorBy]["scale"]
                if isinstance(scale, dict):
                    for value, hex in scale.items():
                        if value not in treeTraits[colorBy]["values"]:
                            warn("Color option \"{}\" specifies a hex code for \"{}\" but this isn't ever seen on the tree nodes.".format(colorBy, value))
                elif isinstance(scale, list):
                    error("List colour scales are invalid")
                else:
                    error("String colour scales are not yet implemented")
            if "domain" in data["colorings"][colorBy]:
                domain = data["colorings"][colorBy]["domain"]
                if data["colorings"][colorBy]["type"] in ["ordinal", "categorical"]:
                    inMetaNotInTree = [val for val in domain if val not in treeTraits[colorBy]["values"]]
                    if len(inMetaNotInTree):
                        warn("Domain for {} defined the following values which are not present on the tree: {}".format(colorBy, inMetaNotInTree.join(", ")))
                    inTreeNotInMeta = [val for val in treeTraits[colorBy]["values"] if val not in domain]
                    if len(inTreeNotInMeta):
                        warn("Tree defined values for {} which were not in the domain: {}".format(colorBy, inTreeNotInMeta.join(", ")))
                elif data["colorings"][colorBy]["type"] == "boolean":
                    error("Cannot povide a domain for a boolean coloring ({})".format(colorBy))
    else:
        warn("No colourings were provided")

    if "filters" in data:
        for filter in data["filters"]:
            if filter not in treeTraits:
                error("The filter \"{}\" does not appear as a property on any tree nodes.".format(filter))

    if "author_info" in data:
        if "authors" not in treeTraits:
            error("\"author_info\" exists in the metadata JSON but \"authors\" are never defined on a tree node attr.")
        else:
            for author in data["author_info"]:
                if author not in treeTraits["authors"]["values"]:
                    warn("Author \"{}\" defined in \"author_info\" but does not appear on any tree node.".format(author))
            for author in treeTraits["authors"]["values"]:
                if author not in data["author_info"]:
                    error("Author \"{}\" defined on a tree node but not in \"author_info\".".format(author))

    genes_with_mutations = collectMutationGenes(data['tree'])
    if len(genes_with_mutations):
        if "genome_annotations" not in data:
            error("The tree defined mutations on genes {}, but annotations aren't defined in the meta JSON.".format(", ".join(genes_with_mutations)))
        else:
            for gene in genes_with_mutations:
                if gene not in data["genome_annotations"]:
                    error("The tree defined mutations on gene {} which doesn't appear in the metadata annotations object.".format(gene))

    return return_status


def register_arguments_v2(subparsers):
    v2 = subparsers.add_parser("v2", help="Validate version 2 JSONs")
    v2.add_argument('--json', required=True, help="JSON to validate")

def run_v2(args):
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
        print("Validating {} using the {} schema (version 2)... ".format(data["path"], schema_type), end='')
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


    if "main" in datasets:
        return_status = verifyMainJSONIsInternallyConsistent(datasets["main"]) | return_status

    if not return_status:
        print("SUCCESS")
    return return_status
