"""
Export JSON files suitable for visualization with Auspice.
The JSON schema is available at <https://nextstrain.org/schemas/dataset/v2>
"""
from .argparse_ import add_command_subparsers
from . import export_v1, export_v2
from textwrap import dedent
import sys

def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("export", help=__doc__)

    # For historical reasons add subparsers for v1 and v2 subcommands:
    subparsers = parser.add_subparsers()
    add_command_subparsers(subparsers, [export_v2, export_v1])

    # Add parsers for the (subcommand-less) "augur export", which is a synonym for "augur export v2"
    export_v2.register_arguments(parser)

    # In augur v28 and earlier, running "augur export" without any extra
    # arguments returned a help-like text instructing us to use a subparser (v1
    # or v2). Now that "augur export" is a synonym for "augur export v2" improve
    # the UX when simply running "augur export" to indicate this change (rather than
    # an error due to the absence of required args)
    error_handler = parser.error
    def error_middleware(message):
        if len(sys.argv)==2 and sys.argv[1]=='export':
            print(dedent("""\
                NOTE: "augur export" is now a synonym for "augur export v2".
                "augur export v1" is deprecated but still valid.
                """), file=sys.stderr)
            parser.print_help()
            sys.exit(2)
        # In all other cases defer to the actual error handler
        error_handler(message)
    parser.error = error_middleware

    return parser

def run(args):
    export_v2.run(args)