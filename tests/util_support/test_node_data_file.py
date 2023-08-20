import json

from augur.__version__ import __version__
from augur.errors import AugurError
from augur.util_support.node_data_file import NodeDataFile

import pytest


@pytest.fixture()
def build_node_data_file(tmpdir):
    def _build_node_data_file(contents):
        with open(f"{tmpdir}/node.json", "w") as file:
            file.write(contents)

        return NodeDataFile(f"{tmpdir}/node.json")

    return _build_node_data_file


class TestNodeDataFile:
    def test_items(self, build_node_data_file):
        node_data_file = build_node_data_file(
            """
            {
                "generated_by": "some app",
                "nodes": { "a": 5 }
            }
            """
        )
        items = node_data_file.items()

        assert len(items) == 1
        assert "generated_by" not in items

    def test_validate_valid(self, build_node_data_file):
        # Implied assertion that no exceptions are raised
        build_node_data_file(
            f"""
            {{
                "annotations": {{ "nuc": {{ "start": 1, "end": 100 }} }},
                "generated_by": {{ "program": "augur", "version": "{__version__}" }},
                "nodes": {{ "a": 5 }}
            }}
            """
        )

    def test_validate_invalid_annotation(self, build_node_data_file):
        with pytest.raises(AugurError, match="annotations.*invalid JSON format"):
            build_node_data_file(
                """
                {
                    "annotations": "hhhhh",
                    "nodes": { "a": 5 }
                }
                """
            )

    def test_validate_nodes_not_dict(self, build_node_data_file):
        with pytest.raises(AugurError, match="is not a dictionary"):
            build_node_data_file('{ "nodes": "hhh" }')

    def test_validate_incompatible_augur_version(self, build_node_data_file):
        with pytest.raises(AugurError, match="incompatibility detected"):
            build_node_data_file(
                """
                {
                    "generated_by": { "program": "augur", "version": "0.1" },
                    "nodes": { "a": 5 }
                }
                """
            )

    def test_validate_malformed_json(self, build_node_data_file):
        with pytest.raises(json.decoder.JSONDecodeError):
            build_node_data_file("hhhhhhhhhh")
