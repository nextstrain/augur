from augur.util_support.node_data import DuplicatedNonDictAttributeError
from augur.util_support.node_data import NodeData

import pytest


class TestNodeData:
    def test_update(self):
        node_data = NodeData()
        node_data.update({"a": 5, "b": {"c": 6}})
        assert node_data.attrs == {"a": 5, "b": {"c": 6}}

        node_data.update({"a": 7})
        assert node_data.attrs == {"a": 7, "b": {"c": 6}}

        node_data.update({"b": {"d": 8}})
        assert node_data.attrs == {"a": 7, "b": {"c": 6, "d": 8}}

        node_data.update({"b": {"d": 9}})
        assert node_data.attrs == {"a": 7, "b": {"c": 6, "d": 9}}

    def test_update_merge_nondict_into_dict(self):
        node_data = NodeData()
        node_data.update({"a": 5, "b": {"c": 6}})

        with pytest.raises(DuplicatedNonDictAttributeError):
            node_data.update({"b": "a string!"})

    def test_update_merge_dict_into_nondict(self):
        node_data = NodeData()
        node_data.update({"a": 5, "b": {"c": 6}})

        with pytest.raises(DuplicatedNonDictAttributeError):
            node_data.update({"a": {"g": 2}})
