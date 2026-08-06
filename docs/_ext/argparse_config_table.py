from docutils import nodes
from sphinx.util.docutils import SphinxDirective
import importlib


def setup(app):
    app.add_directive('argparse-config-table', ArgparseConfigTableDirective)
    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }


class ArgparseConfigTableDirective(SphinxDirective):
    """
    A directive to generate a command-line configuration option table with 2 columns:
    1. Config Option
    2. Description

    Usage:
        .. argparse-config-table:: module.path dict_name
    """
    required_arguments = 2
    optional_arguments = 0
    has_content = False

    def run(self):
        module_path = self.arguments[0]
        dict_name = self.arguments[1]

        module = importlib.import_module(module_path)
        options_dict = getattr(module, dict_name)

        items = []
        for name, spec in options_dict.items():
            desc = spec.get("description", "")
            is_boolean = spec.get("type") == "boolean"

            clean_desc = self._clean_description(desc, is_boolean)
            items.append((name, clean_desc))

        # Create the table structure
        table = nodes.table()
        tgroup = nodes.tgroup(cols=2)
        table.append(tgroup)

        tgroup.append(nodes.colspec(colwidth=30))
        tgroup.append(nodes.colspec(colwidth=70))

        thead = nodes.thead()
        tgroup.append(thead)
        header_row = nodes.row()
        thead.append(header_row)

        header_row.append(nodes.entry('', nodes.paragraph('', 'Config Option')))
        header_row.append(nodes.entry('', nodes.paragraph('', 'Description')))

        tbody = nodes.tbody()
        tgroup.append(tbody)

        for name, desc in items:
            row = nodes.row()
            tbody.append(row)

            opt_entry = nodes.entry()
            opt_entry.append(nodes.literal('', name))
            row.append(opt_entry)

            desc_entry = nodes.entry()
            desc_entry.append(nodes.paragraph('', desc))
            row.append(desc_entry)

        return [table]

    @staticmethod
    def _clean_description(desc: str, is_boolean: bool) -> str:
        desc = desc.replace('\n', ' ').strip()
        desc = desc.replace('%(default).0s', '')
        desc = desc.replace('%(default)s', '')

        if is_boolean:
            desc = f"{desc} (boolean)" if desc else "(boolean)"

        return desc
