import jsonschema
import os
import json
from collections import defaultdict
from pkg_resources import resource_string, resource_filename
from copy import deepcopy

# pip install git+git://github.com/Julian/jsonschema@9632422aa90cb1fbfbbb141954ef6d06437b0801

schemaLocations = {
    "meta": "data/schema_meta.json",
    "tree": "data/schema_tree.json",
    "recommendations": "data/schema_recommendations.json"
}


def loadJSONsToValidate(paths):
    ret = []
    for path in paths:
        with open(path) as f:
            try:
                jsonToValidate = json.load(f)
            except json.JSONDecodeError:
                print("to do JSONDecodeError")
        if path.endswith("_tree.json"):
            name = "tree"
        elif path.endswith("_meta.json"):
            name = "meta"
        else:
            raise Exception('unknown name')
        ret.append({
            "path": path,
            "name": name,
            "json": jsonToValidate
        })
    return ret

def isSchemaValid(schema, name):
    # see http://python-jsonschema.readthedocs.io/en/latest/errors/
    try:
        jsonschema.Draft6Validator.check_schema(schema)
    except jsonschema.exceptions.SchemaError as err:
        print("Schema {} did not pass validation. Error: {}".format(name, err))
        return False
    return True

# with resource_stream(__package__, "data/colors.tsv") as stream:
# resource_string should just return the string
def loadSchemas(names):
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


def checkMetaValuesAppearOnTreeNodes():
    print("checkMetaValuesAppearOnTreeNodes not yet implemented")


def run(args):
    '''
    Validate augur produced JSONs
    '''
    return_code = 0;
    datasets = loadJSONsToValidate(args.json)
    schemas = loadSchemas([x["name"] for x in datasets])

    for data in datasets:
        # make this cmd line controllable
        for schemaType in ["minimal", "recommended"]:
            if schemaType not in schemas[data["name"]]:
                continue
            print("Validating {} using the {}/{} schema... ".format(data["path"], data["name"], schemaType), end='')

            validator = schemas[data["name"]][schemaType]
            try:
                validator.validate(data["json"]) # https://python-jsonschema.readthedocs.io/en/latest/validate/
            except jsonschema.exceptions.ValidationError:
                print("") # new line
                ret_code = 1
                msgType = "ERROR" if schemaType is "minimal" else "WARNING"
                for error in sorted(validator.iter_errors(data["json"]), key=str):
                    # TODO: make this more meaningful. Can you extract the definition if available?
                    print("\t", msgType, error.message, "see:", error.schema_path)
            else:
                print("SUCCESS")

    checkMetaValuesAppearOnTreeNodes()
    return return_code
