class DuplicatedNonDictAttributeError(Exception):
    def __init__(self, key):
        self.key = key


class NodeData:
    """Data structure representing a set of Augur nodes"""

    def __init__(self):
        self.attrs = {}

    def update(self, new_attrs):
        for k, v in new_attrs.items():
            self.deep_add_or_update(self.attrs, k, v)

    def deep_add_or_update(self, d, key, value):
        """
        Recursively merge in new attributes using the following rules:

            if key is new:
                add it
            if key is old, but the old value is not a dict:
                if new value is also not a dict:
                    overwrite it
                if trying to overwrite a non-dict with a dict:
                    raise exception
            if key old and the old value is a dict:
                if new value is also a dict:
                    recursively merge
                if trying to overwrite a dict with a non-dict:
                    raise exception
        """

        if key not in d or (
            not isinstance(d[key], dict) and not isinstance(value, dict)
        ):
            d[key] = value
            return

        if not isinstance(d[key], dict) or not isinstance(value, dict):
            raise DuplicatedNonDictAttributeError(key)

        for sub_k, sub_v in value.items():
            self.deep_add_or_update(d[key], sub_k, sub_v)
