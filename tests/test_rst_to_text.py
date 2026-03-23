"""Tests for RST-to-text conversion used in CLI help output."""
from augur.argparse_ import rst_to_text


class TestRstToText:
    def test_empty(self):
        assert rst_to_text("") == ""
        assert rst_to_text(None) is None

    def test_plain_text(self):
        assert rst_to_text("Hello world.") == "Hello world."

    def test_hyperlink(self):
        text = "Inferred using `TreeTime <https://example.com>`_."
        assert rst_to_text(text) == "Inferred using TreeTime(https://example.com)."

    def test_bold_to_uppercase(self):
        assert rst_to_text("**Comparison methods**") == "COMPARISON METHODS"

    def test_note_directive(self):
        text = """\
        Some text.

        .. note::

            The mutation positions are one-based.
        """
        result = rst_to_text(text)
        assert "Note: The mutation positions are one-based." in result
        assert ".. note::" not in result

    def test_code_block_directive(self):
        text = """\
        Example:

        .. code-block:: json

            {
                "default": 1,
                "map": {}
            }
        """
        result = rst_to_text(text)
        assert ".. code-block::" not in result
        assert '"default": 1' in result
        assert '"map": {}' in result

    def test_numbered_list(self):
        text = """\
        Methods include:

        #. root: all nodes
        #. ancestor: each tip
        #. pairwise: all tips
        """
        result = rst_to_text(text)
        assert "1. root: all nodes" in result
        assert "2. ancestor: each tip" in result
        assert "3. pairwise: all tips" in result
        assert "#." not in result

    def test_single_backtick_removal(self):
        assert rst_to_text("`augur translate`") == "augur translate"

    def test_double_backtick_preserved(self):
        assert rst_to_text("``literal``") == "``literal``"

    def test_combined_rst(self):
        """Test a docstring with multiple RST constructs (like ancestral.py)."""
        text = """\
        Infer ancestral sequences based on a tree.

        The ancestral sequences are inferred using `TreeTime <https://example.com>`_.
        Each node gets assigned a list of mutations.

        Use `augur translate` for amino acid mutations.

        .. note::

            The mutation positions are one-based.
        """
        result = rst_to_text(text)
        assert "TreeTime(https://example.com)" in result
        assert "augur translate" in result
        assert "`" not in result
        assert ".. note::" not in result
        assert "Note: The mutation positions are one-based." in result
