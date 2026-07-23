from docutils import nodes
from sphinx.util.docutils import SphinxDirective
import importlib
import argparse

def setup(app):
    app.add_directive('argparse-config-table', ArgparseConfigTableDirective)
    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }

class ArgparseConfigTableDirective(SphinxDirective):
    """
    A directive to generate a table with 2 columns:
    1. Config Option
    2. Description

    Usage:
        .. argparse-config-table:: module.path function_name
    """
    required_arguments = 2
    optional_arguments = 0
    has_content = False

    def run(self):
        module_path = self.arguments[0]
        func_name = self.arguments[1]

        module = importlib.import_module(module_path)
        func = getattr(module, func_name)
        
        # Instantiate the parser (assumes it takes a subparsers object)
        p = argparse.ArgumentParser()
        sp = p.add_subparsers()
        parser = func(sp)

        # Create the table structure
        table = nodes.table()
        tgroup = nodes.tgroup(cols=2)
        table.append(tgroup)

        # Define column specifications
        tgroup.append(nodes.colspec(colwidth=30))
        tgroup.append(nodes.colspec(colwidth=70))

        # Create table header
        thead = nodes.thead()
        tgroup.append(thead)
        header_row = nodes.row()
        thead.append(header_row)

        header_row.append(nodes.entry('', nodes.paragraph('', 'Config Option')))
        header_row.append(nodes.entry('', nodes.paragraph('', 'Description')))

        # Create table body
        tbody = nodes.tbody()
        tgroup.append(tbody)

        for action in parser._actions:
            if action.dest in ('help', 'config'):
                continue
            
            name = action.dest
            
            # Clean up the help string
            desc = action.help if action.help else ''
            desc = desc.replace('\n', ' ').strip()
            desc = desc.replace('%(default).0s', '')
            desc = desc.replace('%(default)s', '')

            # Indicate boolean options
            if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
                if desc:
                    desc += " (boolean)"
                else:
                    desc = "(boolean)"

            row = nodes.row()
            tbody.append(row)

            opt_entry = nodes.entry()
            opt_entry.append(nodes.literal('', name))
            row.append(opt_entry)

            desc_entry = nodes.entry()
            desc_entry.append(nodes.paragraph('', desc))
            row.append(desc_entry)

        return [table]
