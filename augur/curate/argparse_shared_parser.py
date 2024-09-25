import argparse
from augur.argparse_ import ExtendOverwriteDefault
from augur.io.metadata import DEFAULT_DELIMITERS
from augur.types import DataErrorMethod


def create_shared_parser():
    """
    Creates an argparse.ArgumentParser that is intended to be used as a parent
    parser¹ for all `augur curate` subcommands. This should include all options
    that are intended to be shared across the subcommands.

    Note that any options strings used here cannot be used in individual subcommand
    subparsers unless the subparser specifically sets `conflict_handler='resolve'` ²,
    then the subparser option will override the option defined here.

    Based on https://stackoverflow.com/questions/23296695/permit-argparse-global-flags-after-subcommand/23296874#23296874

    ¹ https://docs.python.org/3/library/argparse.html#parents
    ² https://docs.python.org/3/library/argparse.html#conflict-handler
    """
    shared_parser = argparse.ArgumentParser(add_help=False)

    shared_inputs = shared_parser.add_argument_group(
        title="INPUTS",
        description="""
            Input options shared by all `augur curate` commands.
            If no input options are provided, commands will try to read NDJSON records from stdin.
        """)
    shared_inputs.add_argument("--metadata",
        help="Input metadata file. May be plain text (TSV, CSV) or an Excel or OpenOffice spreadsheet workbook file. When an Excel or OpenOffice workbook, only the first visible worksheet will be read and initial empty rows/columns will be ignored. Accepts '-' to read plain text from stdin.")
    shared_inputs.add_argument("--id-column",
        help="Name of the metadata column that contains the record identifier for reporting duplicate records. "
             "Uses the first column of the metadata file if not provided. "
             "Ignored if also providing a FASTA file input.")
    shared_inputs.add_argument("--metadata-delimiters", default=DEFAULT_DELIMITERS, nargs="+", action=ExtendOverwriteDefault,
        help="Delimiters to accept when reading a plain text metadata file. Only one delimiter will be inferred.")

    shared_inputs.add_argument("--fasta",
        help="Plain or gzipped FASTA file. Headers can only contain the sequence id used to match a metadata record. " +
             "Note that an index file will be generated for the FASTA file as <filename>.fasta.fxi")
    shared_inputs.add_argument("--seq-id-column",
        help="Name of metadata column that contains the sequence id to match sequences in the FASTA file.")
    shared_inputs.add_argument("--seq-field",
        help="The name to use for the sequence field when joining sequences from a FASTA file.")

    shared_inputs.add_argument("--unmatched-reporting",
        type=DataErrorMethod.argtype,
        choices=list(DataErrorMethod),
        default=DataErrorMethod.ERROR_FIRST,
        help="How unmatched records from combined metadata/FASTA input should be reported.")
    shared_inputs.add_argument("--duplicate-reporting",
        type=DataErrorMethod.argtype,
        choices=list(DataErrorMethod),
        default=DataErrorMethod.ERROR_FIRST,
        help="How should duplicate records be reported.")

    shared_outputs = shared_parser.add_argument_group(
        title="OUTPUTS",
        description="""
            Output options shared by all `augur curate` commands.
            If no output options are provided, commands will output NDJSON records to stdout.
        """)
    shared_outputs.add_argument("--output-metadata",
        help="Output metadata TSV file. Accepts '-' to output TSV to stdout.")

    shared_outputs.add_argument("--output-fasta",
        help="Output FASTA file.")
    shared_outputs.add_argument("--output-id-field",
        help="The record field to use as the sequence identifier in the FASTA output.")
    shared_outputs.add_argument("--output-seq-field",
        help="The record field that contains the sequence for the FASTA output. "
             "This field will be deleted from the metadata output.")

    return shared_parser


shared_parser = create_shared_parser()
