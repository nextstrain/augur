import jsonschema
import os, sys
import json
from collections import defaultdict
from pkg_resources import resource_string, resource_filename
from copy import deepcopy

schemaLocations = {
    "meta": "data/schema_meta.json",
    "tree": "data/schema_tree.json",
    "recommendations": "data/schema_recommendations.json"
}

def loadJSONsToValidate(paths):
    """
    paths: array of paths to JSONs. Must end in "_tree.json", "_meta.json" etc.
    """
    ret = {}
    for path in paths:
        with open(path) as f:
            try:
                jsonToValidate = json.load(f)
            except json.JSONDecodeError:
                raise(Exception("Input JSON Decode Error"))
        if path.endswith("_tree.json"):
            name = "tree"
        elif path.endswith("_meta.json"):
            name = "meta"
        else:
            raise Exception('Unknown JSON type for {}'.format(path))
        if name in ret:
            raise Exception('Cannot deal with multiple JSONs of the same type ({})'.format(name))
        ret[name] = {
            "path": path,
            "json": jsonToValidate
        }
    return ret

def isSchemaValid(schema, name):
    # see http://python-jsonschema.readthedocs.io/en/latest/errors/
    try:
        jsonschema.Draft6Validator.check_schema(schema)
    except jsonschema.exceptions.SchemaError as err:
        print("Schema {} did not pass validation. Error: {}".format(name, err))
        return False
    return True

def loadSchemas(names):
    """
    For names (such as "tree", "meta"), load and internally-verify the schema.
    If there are additional fields that are recommended (see the file specified by schemaLocations["recommendations"])
    then dynamically modify the schema to create a (second) recommended schema.
    """
    ret = defaultdict(dict)

    # Unfortunately, you cannot create a recommended schemas which inherit the minimal schema using $ref
    # see https://spacetelescope.github.io/understanding-json-schema/reference/combining.html
    # so we do it dynamically here
    try:
        recommendations = json.loads(resource_string(__package__, schemaLocations["recommendations"]))
    except:
        print("Could not load the schema recommendations. Proceeding with minimal schemas only.")

    for name in names:
        if name in ret:
            continue
        if name not in schemaLocations:
            print("ERROR: No specified schema for ", name)
            continue

        ### LOAD THE (MINIMAL) SCHEMA
        try:
            schema = json.loads(resource_string(__package__, schemaLocations[name]))
        except json.JSONDecodeError as err:
            print("Schema {} is not a valid JSON file. Error: {}".format(schemaLocations[name], err))
            continue
        if not isSchemaValid(schema, schemaLocations[name]):
            continue
        ret[name]["minimal"] = jsonschema.Draft6Validator(schema)

        # dynamically create the "recommended" schema
        if recommendations and name in recommendations:
            schema = deepcopy(schema) # so we don't modify the minimal one
            for keyToInject, value in recommendations[name].items():
                schema[keyToInject] = value;
            if not isSchemaValid(schema, "Recommended schema derived from "+schemaLocations[name]):
                continue
            ret[name]["recommended"] = jsonschema.Draft6Validator(schema)

    return ret

def collectTreeAttrs(root):
    """
    Collect all keys specified on node->attr throughout the tree
    If the values of these keys are strings, then also collect the values
    """
    seen = defaultdict(lambda: {"count": 0, "values": set(), "onAllNodes": False})
    num_nodes, num_terminal = (0, 0)
    def recurse(node):
        nonlocal num_nodes, num_terminal
        num_nodes += 1
        for property in node["attr"].keys():
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


def verifyMetaJSONIsInternallyConsistent(status, meta):
    """
    Check all possible sources of internal conflict within the metadata JSON that cannot be checked by the schema
    """
    print("Validating that {} is internally consistent... ".format(meta["path"]))
    j = meta["json"]
    if "panels" in j and "entropy" in j["panels"] and "annotations" not in j:
        status = "error"
        print("\tERROR: The entropy panel has been specified but annotations don't exist.", file=sys.stderr)

    return(status)


def verifyMetaAndTreeJSONsAreInternallyConsistent(status, meta, tree):
    """
    Check all possible sources of conflict between the metadata & tree JSONs that cannot be checked by the schemas
    """
    print("Validating that {} and {} are internally consistent... ".format(meta["path"], tree["path"]))
    mj, tj = (meta["json"], tree["json"])
    treeAttrs, num_terminal_nodes = collectTreeAttrs(tj)

    def error(msg):
        nonlocal status
        status = "error"
        print("\tERROR: ", msg, file=sys.stderr)

    def warn(msg):
        nonlocal status
        if status == "ok":
            status = "warn"
        print("\tWARNING: ", msg, file=sys.stderr)


    if "geo" in mj:
        for geoName in mj["geo"]:
            if geoName not in treeAttrs:
                error("The geographic resolution \"{}\" does not appear as an attr on any tree nodes.".format(geoName))
            else:
                for geoValue in mj["geo"][geoName]:
                    if geoValue not in treeAttrs[geoName]["values"]:
                        warn("\"{}\", a value of the geographic resolution \"{}\", does not appear as a value of attr->{} on any tree nodes.".format(geoValue, geoName, geoName))

    if "color_options" in mj:
        for colorBy in [x for x in mj["color_options"] if x != "gt"]:
            if colorBy not in treeAttrs:
                error("The color_option \"{}\" does not appear as an attr on any tree nodes.".format(colorBy))
            elif "color_map" in mj["color_options"][colorBy]:
                for (value, hex) in mj["color_options"][colorBy]["color_map"]:
                    if value not in treeAttrs[colorBy]["values"]:
                        warn("Color option \"{}\" specifies a hex code for \"{}\" but this isn't ever seen on the tree nodes.".format(colorBy, value))

    if "filters" in mj:
        for filter in mj["filters"]:
            if filter not in treeAttrs:
                error("The filter \"{}\" does not appear as an attr on any tree nodes.".format(filter))

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

    return status


def run(args):
    '''
    Validate auspice-compatable JSONs against a schema
    '''
    status = "ok" # possible values are "warn" and "error"

    try:
        datasets = loadJSONsToValidate(args.json)
        schemas = loadSchemas([x for x in datasets])
    except Exception as err:
        print(err, file=sys.stderr)
        return 1

    for name, data in datasets.items():
        # make this cmd line controllable
        for schemaType in ["minimal", "recommended"]:
            if schemaType not in schemas[name]:
                continue
            print("Validating {} using the {}/{} schema... ".format(data["path"], name, schemaType), end='')

            validator = schemas[name][schemaType]
            try:
                validator.validate(data["json"]) # https://python-jsonschema.readthedocs.io/en/latest/validate/
            except jsonschema.exceptions.ValidationError:
                print("FAILED")
                if schemaType == "minimal":
                    status = "error"
                    msgType = "ERROR"
                elif status == "ok":
                    status = "warn"
                    msgType = "WARNING"

                for error in sorted(validator.iter_errors(data["json"]), key=str):
                    trace = list(error.schema_path) # this may be really long for nested tree structures
                    if len(trace) > 6:
                        trace = ["..."] + trace[-5:]
                    trace = [str(x) for x in trace]
                    print("\t{}: {}. Trace: {}".format(msgType, error.message, " - ".join(trace)), file=sys.stderr)
            else:
                print("SUCCESS")

    if "meta" in datasets:
        status = verifyMetaJSONIsInternallyConsistent(status, datasets["meta"])
    if "meta" in datasets and "tree" in datasets:
        status = verifyMetaAndTreeJSONsAreInternallyConsistent(status, datasets["meta"], datasets["tree"])


    return 0 if status != "error" else 1
