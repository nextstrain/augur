"""
Custom Sphinx directives for augur subsample documentation.
"""
from docutils import nodes
from sphinx.util.docutils import SphinxDirective
from docutils.statemachine import StringList
import importlib
import json
import os


def setup(app):
    app.add_directive('yaml-option-table', YAMLOptionTableDirective)
    app.add_directive('cli-option-table', CLIOptionTableDirective)
    app.add_directive('sample-options-schema-table', SampleOptionsSchemaTableDirective)

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

        # Import the module and get the configuration dictionary
        module = importlib.import_module(module_path)
        config_dict = getattr(module, var_name)

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
        # Merge all configuration dictionaries
        config_dict = {}
        for dict_path in self.arguments:
            module_path, var_name = dict_path.rsplit('.', 1)

            # Import the module and get the configuration dictionary
            module = importlib.import_module(module_path)
            single_config_dict = getattr(module, var_name)
            config_dict.update(single_config_dict)

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

        for key in config_dict.keys():
            value = config_dict[key]
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


class SampleOptionsSchemaTableDirective(SphinxDirective):
    """
    A directive to generate a table with sample options from the JSON schema.

    Creates a table with columns:
    1. Property - YAML key name
    2. Type - Data type information
    3. Description - Human-readable description

    Usage:
        .. sample-options-schema-table:: path/to/schema.json
    """
    required_arguments = 1
    optional_arguments = 0
    has_content = False

    def run(self):
        schema_path = self.arguments[0]

        # Open the schema file relative to the augur root
        base_path = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        full_schema_path = os.path.join(base_path, schema_path)
        with open(full_schema_path, 'r') as f:
            schema = json.load(f)

        # Extract sample options from schema
        sample_options = schema['$defs']['sampleProperties']['properties']

        # Create the table structure
        table = nodes.table()
        tgroup = nodes.tgroup(cols=3)
        table.append(tgroup)

        # Define column specifications
        tgroup.append(nodes.colspec(colwidth=1))
        tgroup.append(nodes.colspec(colwidth=1))
        tgroup.append(nodes.colspec(colwidth=2))

        # Create table header
        thead = nodes.thead()
        tgroup.append(thead)
        header_row = nodes.row()
        thead.append(header_row)

        header_row.append(nodes.entry('', nodes.paragraph('', 'Option')))
        header_row.append(nodes.entry('', nodes.paragraph('', 'Type')))
        header_row.append(nodes.entry('', nodes.paragraph('', 'Description')))

        # Create table body
        tbody = nodes.tbody()
        tgroup.append(tbody)

        for prop_name in sample_options.keys():
            prop_def = sample_options[prop_name]
            row = nodes.row()
            tbody.append(row)

            # Option column
            option_entry = nodes.entry()
            option_entry.append(nodes.literal('', prop_name))
            row.append(option_entry)

            # Type column
            type_entry = nodes.entry()
            type_info = self._format_type_info(prop_def)
            type_entry.append(type_info)
            row.append(type_entry)

            # Description column
            description_entry = nodes.entry()
            description_entry['classes'].append('wrap-children-60ch')
            description = prop_def.get('description', '')
            container = nodes.container()
            self.state.nested_parse(StringList(description.split('\n')), 0, container)
            for child in container.children:
                description_entry += child
            row.append(description_entry)

        return [table]

    def _format_type_info(self, prop_def):
        """Format type information from JSON schema property definition."""
        prop_type = prop_def.get('type', 'unknown')

        if prop_type == 'array':
            items_type = prop_def.get('items', {}).get('type', 'unknown')
            return nodes.paragraph('', f"list of {items_type}s")
        elif prop_type == 'string' and 'enum' in prop_def:
            container = nodes.container()
            container += nodes.paragraph('', 'one of:')
            bullet_list = nodes.bullet_list()
            for value in prop_def['enum']:
                list_item = nodes.list_item()
                list_item += nodes.literal('', value)
                bullet_list += list_item
            container += bullet_list
            return container
        else:
            return nodes.paragraph('', prop_type)
