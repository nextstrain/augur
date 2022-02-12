import functools

from augur.validation.json_schema import JsonSchema


class JsonValidator:
    def __init__(self, json_obj, schema_name):
        self.json_obj = json_obj
        self.schema_name = schema_name

    @property
    def schema(self):
        return JsonSchema(self.schema_name)

    @property
    @functools.lru_cache()
    def errors(self):
        return self.schema.validation_errors(self.json_obj)

    @property
    def is_valid(self):
        return len(self.errors) == 0
