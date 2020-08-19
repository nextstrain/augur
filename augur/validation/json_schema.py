import functools
import json
import jsonschema
from pkg_resources import resource_string


def truncate_trace(trace):
    if len(trace) <= 6:
        return trace

    return ["..."] + trace[-5:]


class JsonSchema:
    def __init__(self, name):
        self.name = name

        self.validate_schema()

    def validation_errors(self, json):
        """Validate provided JSON against the schema"""
        return [
            f"\tERROR: {error.message}. Trace: {' - '.join(truncate_trace(error.schema_path))}"
            for error in sorted(self.schema.iter_errors(json), key=str)
        ]

    @property
    @functools.lru_cache()
    def schema(self):
        return jsonschema.Draft6Validator(self.schema_src)

    @property
    @functools.lru_cache()
    def schema_src(self):
        return json.loads(resource_string("augur.data", self.name))

    def validate_schema(self):
        """Validate the schema itself"""
        jsonschema.Draft6Validator.check_schema(self.schema_src)
