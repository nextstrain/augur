import jsonschema
import os
import json
from collections import defaultdict
from pkg_resources import resource_string

# pip install git+git://github.com/Julian/jsonschema@9632422aa90cb1fbfbbb141954ef6d06437b0801


schemas = {
    "meta": {
        "minimal": "data/schema_meta_minimal.json",
        "recommended": None
    },
    "tree": {
        "minimal": "data/schema_tree_minimal.json",
        "recommended": None
    },
    "frequencies": {
        "minimal": None,
        "recommended": None
    }
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

# with resource_stream(__package__, "data/colors.tsv") as stream:
# resource_string should just return the string
def loadSchemas(names):
    ret = defaultdict(dict)
    for name in names:
        if name in ret:
            continue
        for schemaType, schemaPath in schemas[name].items():
            if schemaPath is None:
                continue

            # see http://python-jsonschema.readthedocs.io/en/latest/errors/
            try:
                schema = json.loads(resource_string(__package__, schemaPath))
                jsonschema.Draft6Validator.check_schema(schema)
            except json.JSONDecodeError as err:
                print("Schema {} is not a valid JSON file. Error: {}".format(schemaPath, err))
            except jsonschema.exceptions.SchemaError as err:
                print("Schema {} did not pass validation. Error: {}".format(schemaPath, err))
            ret[name][schemaType] = schema

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
                # print("WARNING: {} schema for {} doesn't exist!".format(schemaType, data["name"]))
                continue;
            print("Validating {} using the {}/{} schema".format(data["path"], data["name"], schemaType))

            schema = schemas[data["name"]][schemaType]

            # https://python-jsonschema.readthedocs.io/en/latest/validate/
            try:
                jsonschema.Draft6Validator(schema).validate(data["json"])
            except jsonschema.exceptions.ValidationError:
                ret_code = 1
                v = jsonschema.Draft6Validator(schema)
                for error in sorted(v.iter_errors(data["json"]), key=str):
                    print("\tERROR", error.message, "see:", error.schema_path)
            else:
                print("\tSuccess")

    checkMetaValuesAppearOnTreeNodes()
    return return_code
