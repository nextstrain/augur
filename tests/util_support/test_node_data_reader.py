from augur.util_support.node_data_reader import NodeDataReader

import pytest


@pytest.fixture()
def prepare_file(tmpdir):
    def _prepare_file(name, contents):
        with open(f"{tmpdir}/{name}", "w") as file:
            file.write(contents)

    return _prepare_file


class TestNodeDataFile:
    def test_read(self, tmpdir, prepare_file):
        prepare_file(
            "file1.json",
            """
            {
                "a": 5,
                "nodes": {
                    "NODE_1": { "x": 5 },
                    "NODE_2": { "x": 6 },
                    "NODE_3": { "x": 7 }
                }
            }
            """,
        )
        prepare_file(
            "file2.json",
            """
            {
                "nodes": {
                    "NODE_3": { "y": 8 },
                    "NODE_4": { "x": 9 }
                }
            }
            """,
        )

        node_data = NodeDataReader(
            [f"{tmpdir}/file1.json", f"{tmpdir}/file2.json"]
        ).read()
        assert node_data == {
            "a": 5,
            "nodes": {
                "NODE_1": {"x": 5},
                "NODE_2": {"x": 6},
                "NODE_3": {"x": 7, "y": 8},
                "NODE_4": {"x": 9},
            },
        }

    def test_read_bad_file(self):
        NodeDataReader(["/does/not/exist.json"])

    def test_read_dict_nonuniformity(self, tmpdir, prepare_file):
        prepare_file(
            "file1.json",
            """
            {
                "nodes": {"node_name": "some_value"},
                "a": {}
            }
            """,
        )
        prepare_file(
            "file2.json",
            """
            {
                "nodes": {"node_name": "some_other_value"},
                "a": "nah"
            }
            """,
        )
        with pytest.raises(SystemExit):
            NodeDataReader([f"{tmpdir}/file1.json", f"{tmpdir}/file2.json"]).read()

    def test_read_check_against_tree(self, tmpdir, prepare_file):
        prepare_file(
            "file1.json",
            """
            {
                "a": 5,
                "nodes": {
                    "NODE_1": { "x": 5 },
                    "NODE_2": { "x": 6 },
                    "NODE_3": { "x": 7 }
                }
            }
            """,
        )
        prepare_file("tree.newick", "(NODE_1, NODE_2) NODE_3")

        NodeDataReader([f"{tmpdir}/file1.json"], f"{tmpdir}/tree.newick").read()

    def test_read_check_against_tree_bad(self, tmpdir, prepare_file):
        prepare_file(
            "file1.json",
            """
            {
                "a": 5,
                "nodes": {
                    "NODE_1": { "x": 5 },
                    "NODE_2": { "x": 6 },
                    "NODE_3": { "x": 7 }
                }
            }
            """,
        )
        prepare_file("tree.newick", "Noooooooooope")

        with pytest.raises(SystemExit):
            NodeDataReader([f"{tmpdir}/file1.json"], f"{tmpdir}/tree.newick").read()

    def test_read_check_against_missing_tree(self, tmpdir):
        with pytest.raises(SystemExit):
            node_names_from_tree = NodeDataReader(
                [f"{tmpdir}/file1.json"],
                f"{tmpdir}/missing_file.txt"
            ).node_names_from_tree_file
