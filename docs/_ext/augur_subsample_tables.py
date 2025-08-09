"""
Custom Sphinx directives for augur subsample documentation.
"""
from docutils import nodes
from sphinx.util.docutils import SphinxDirective
import importlib


def setup(app):
    app.add_directive('yaml-option-table', YAMLOptionTableDirective)
    app.add_directive('cli-option-table', CLIOptionTableDirective)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }


class YAMLOptionTableDirective(SphinxDirective):
    """
    A directive to generate a table with 2 columns:

    1. YAML config option
    2. augur filter CLI option

    Usage:
        .. yaml-option-table:: module.path.DICT_NAME
    """
    required_arguments = 1
    optional_arguments = 0
    has_content = False

    def run(self):
        # Parse the module path and variable name
        module_path, var_name = self.arguments[0].rsplit('.', 1)

        try:
            # Import the module and get the configuration dictionary
            module = importlib.import_module(module_path)
            config_dict = getattr(module, var_name)
        except (ImportError, AttributeError) as e:
            error_msg = f"Could not import {self.arguments[0]}: {e}"
            error_node = nodes.error('', nodes.paragraph('', error_msg))
            return [error_node]

        # Create the table structure
        table = nodes.table()
        tgroup = nodes.tgroup(cols=2)
        table.append(tgroup)

        # Define column specifications
        tgroup.append(nodes.colspec(colwidth=1))
        tgroup.append(nodes.colspec(colwidth=1))

        # Create table header
        thead = nodes.thead()
        tgroup.append(thead)
        header_row = nodes.row()
        thead.append(header_row)

        header_row.append(nodes.entry('', nodes.paragraph('', 'YAML config option')))
        header_row.append(nodes.entry('', nodes.paragraph('', 'augur filter CLI option')))

        # Create table body
        tbody = nodes.tbody()
        tgroup.append(tbody)

        for key in config_dict.keys():
            value = config_dict[key]
            row = nodes.row()
            tbody.append(row)

            yaml_option = nodes.entry()
            yaml_option.append(nodes.literal('', key))
            row.append(yaml_option)

            filter_option = nodes.entry()
            if isinstance(value, tuple):
                if value[1] is None:
                    filter_option.append(nodes.literal('', value[0]))
                else:
                    filter_option.append(nodes.literal('', value[0]))
                    filter_option.append(nodes.Text(' / '))
                    filter_option.append(nodes.literal('', value[1]))
            else:
                # Simple string flag
                filter_option.append(nodes.literal('', value))

            row.append(filter_option)

        return [table]


class CLIOptionTableDirective(SphinxDirective):
    """
    A directive to generate a table with 2 columns:

    1. augur subsample CLI option
    2. augur filter CLI option

    Usage:
        .. cli-option-table:: module.path.DICT_NAME [module.path.DICT_NAME2 ...]
    """
    required_arguments = 1
    optional_arguments = 10  # Allow multiple dictionary arguments
    has_content = False

    def run(self):
        # Collect all dictionaries to merge
        all_dicts = {}

        # Process all provided dictionary arguments
        for dict_path in self.arguments:
            module_path, var_name = dict_path.rsplit('.', 1)

            try:
                # Import the module and get the configuration dictionary
                module = importlib.import_module(module_path)
                config_dict = getattr(module, var_name)
                all_dicts.update(config_dict)
            except (ImportError, AttributeError) as e:
                error_msg = f"Could not import {dict_path}: {e}"
                error_node = nodes.error('', nodes.paragraph('', error_msg))
                return [error_node]

        # Create the table structure
        table = nodes.table()
        tgroup = nodes.tgroup(cols=2)
        table.append(tgroup)

        # Define column specifications
        tgroup.append(nodes.colspec(colwidth=1))
        tgroup.append(nodes.colspec(colwidth=1))

        # Create table header
        thead = nodes.thead()
        tgroup.append(thead)
        header_row = nodes.row()
        thead.append(header_row)

        header_row.append(nodes.entry('', nodes.paragraph('', 'augur subsample CLI option')))
        header_row.append(nodes.entry('', nodes.paragraph('', 'augur filter CLI option')))

        # Create table body
        tbody = nodes.tbody()
        tgroup.append(tbody)

        for key in all_dicts.keys():
            value = all_dicts[key]
            row = nodes.row()
            tbody.append(row)

            subsample_option = nodes.entry()
            subsample_flag = f"--{key.replace('_', '-')}"
            subsample_option.append(nodes.literal('', subsample_flag))
            row.append(subsample_option)

            filter_option = nodes.entry()
            if isinstance(value, tuple):
                if len(value) == 2 and value[1] is None:
                    filter_option.append(nodes.literal('', value[0]))
                elif len(value) == 2:
                    filter_option.append(nodes.literal('', value[0]))
                    filter_option.append(nodes.Text(' / '))
                    filter_option.append(nodes.literal('', value[1]))
                else:
                    filter_option.append(nodes.literal('', str(value)))
            else:
                filter_option.append(nodes.literal('', value))
            row.append(filter_option)

        return [table]
