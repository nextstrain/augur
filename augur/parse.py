"""
Parse delimited fields from FASTA sequence names into a TSV and FASTA file.
"""
import Bio.SeqRecord
import pandas as pd
import sys
from typing import Dict, Sequence, Tuple

from .argparse_ import ExtendOverwriteDefault
from .io.file import open_file
from .io.sequences import read_sequences, write_sequences
from .dates import get_numerical_date_from_value
from .errors import AugurError

PARSE_DEFAULT_ID_COLUMNS = ("strain", "name")

forbidden_characters = str.maketrans(
    {' ': None,
     '(': '_',
     ')': '_',
     '[': '_',
     ']': '_',
     ':': '_',
     ',': '_',
     ';': '_',
     '\\': '_'}
)

def fix_dates(d: str, dayfirst: bool = True) -> str:
    '''
    attempt to parse a date string using pandas date parser. If ambiguous,
    the argument 'dayfirst' determines whether month or day is assumed to be
    the first field. Incomplete dates will be padded with XX.
    On failure to parse the date, the function will return the input.
    '''
    try:
        try:
            # pandas <= 2.1
            from pandas.core.tools.datetimes import parsing  # type: ignore[attr-defined, import-not-found]
        except ImportError:
            # pandas >= 2.2
            from pandas._libs.tslibs import parsing
        try:
            # pandas 2.x
            results = parsing.parse_datetime_string_with_reso(d, dayfirst=dayfirst)
        except AttributeError:
            # pandas 1.x
            results = parsing.parse_time_string(d, dayfirst=dayfirst)  # type: ignore[attr-defined]
        if len(results) == 2:
            dto, res = results
        else:
            dto, _, res = results

        if res == 'year':
            return "%d-XX-XX"%dto.year
        elif res == 'month':
            return "%d-%02d-XX"%(dto.year, dto.month)
        else:
            return "%d-%02d-%02d"%(dto.year, dto.month, dto.day)
    except ValueError as e:
        # If the date can't be parsed by pandas above or as our own ambiguous
        # date format (e.g., "2020-XX-XX"), let the user know.
        try:
            parsed_date = get_numerical_date_from_value(d, "%Y-%m-%d")
        except ValueError:
            parsed_date = None

        if parsed_date is None:
            print("WARNING: unable to parse %s as date"%d, e, file=sys.stderr)

        return d

def prettify(x, trim=0, camelCase=False, etal=None, removeComma=False):
    res = x
    if 0 < trim < len(x):
        res = x[:trim] + "..."

    if res in {'usvi', 'usa', 'uk'}:
        res = res.upper()

    words = res.split('_')

    if camelCase:
        words = map(str.capitalize, words)

    res = ' '.join(words)

    if removeComma:
        res = res.replace(',', '')

    if etal == 'lower':
        res = res.replace('Et Al', 'et al').replace('Et Al.', 'et al.').replace('Et al', 'et al')
    elif etal == 'strip':
        res = res.replace('et al.', '').replace('Et Al.', '').replace('et al', '').replace('Et Al', '')

    return res


def parse_sequence(
        sequence: Bio.SeqRecord.SeqRecord,
        fields: Sequence[str],
        strain_key: str,
        separator: str,
        prettify_fields: Sequence[str],
        fix_dates_format: str,
    ) -> Tuple[Bio.SeqRecord.SeqRecord, Dict[str, str]]:
    """Parse a single sequence record into a sequence record and associated metadata.

    Parameters
    ----------
    sequence
        a BioPython sequence record to parse with metadata stored in its description field.

    fields
        a list of names for fields expected in the given record's description.

    strain_key
        name of the field to use as the given sequence's unique id

    separator
        delimiter to split record description by.

    prettify_fields
        a list of field names for which the values in those fields should be prettified.

    fix_dates_format
        parse "date" field into the requested canonical format ("dayfirst" or "monthfirst").

    Returns
    -------
    Bio.SeqRecord.SeqRecord :
        a BioPython sequence record with the given sequence's name as the record
        id and all other metadata stripped.

    dict :
        metadata associated with the given record indexed by the given field names.
    """
    sequence_fields = map(str.strip, sequence.description.split(separator))
    metadata = dict(zip(fields, sequence_fields))

    metadata[strain_key] = metadata[strain_key].translate(forbidden_characters)
    sequence.name = sequence.id = metadata[strain_key]
    sequence.description = ''

    if prettify_fields:
        for field in metadata.keys() & prettify_fields:
            if isinstance(metadata[field], str):
                metadata[field] = prettify(metadata[field], camelCase=(not field.startswith('author')),
                                            etal='lower' if field.startswith('author') else None)

    # parse dates and convert to a canonical format
    if fix_dates_format and 'date' in metadata:
        metadata['date'] = fix_dates(
            metadata['date'],
            dayfirst=fix_dates_format=='dayfirst'
        )

    return sequence, metadata


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("parse", help=__doc__)
    parser.add_argument('--sequences', '-s', required=True, help="sequences in fasta or VCF format")
    parser.add_argument('--output-sequences', required=True, help="output sequences file")
    parser.add_argument('--output-metadata', required=True, help="output metadata file")
    parser.add_argument('--output-id-field', required=False,
                        help=f"The record field to use as the sequence identifier in the FASTA output. If not provided, this will use the first available of {PARSE_DEFAULT_ID_COLUMNS}. If none of those are available, this will use the first field in the fasta header.")
    parser.add_argument('--fields', required=True, nargs='+', action=ExtendOverwriteDefault, help="fields in fasta header")
    parser.add_argument('--prettify-fields', nargs='+', action=ExtendOverwriteDefault, help="apply string prettifying operations (underscores to spaces, capitalization, etc) to specified metadata fields")
    parser.add_argument('--separator', default='|', help="separator of fasta header")
    parser.add_argument('--fix-dates', choices=['dayfirst', 'monthfirst'],
                                help="attempt to parse non-standard dates and output them in standard YYYY-MM-DD format")
    return parser


def run(args):
    '''
    parse a fasta file and turn information in the header into
    a tsv or csv file.
    '''
    sequences = read_sequences(args.sequences)

    # if strain or name are found in specified fields, use this
    # field to index the dictionary and the data frame
    meta_data = {}

    strain_key = None
    if args.output_id_field:
        if args.output_id_field not in args.fields:
            raise AugurError(f"Output id field '{args.output_id_field}' not found in fields {args.fields}.")
        strain_key = args.output_id_field
    else:
        for possible_id in PARSE_DEFAULT_ID_COLUMNS:
            if possible_id in args.fields:
                strain_key = possible_id
                break
        if not strain_key:
            strain_key = args.fields[0]

    # loop over sequences, parse fasta header of each sequence
    with open_file(args.output_sequences, "wt") as handle:
        for sequence in sequences:
            sequence_record, sequence_metadata = parse_sequence(
                sequence,
                args.fields,
                strain_key,
                args.separator,
                args.prettify_fields,
                args.fix_dates
            )
            if sequence_record.id in meta_data:
                raise AugurError(f"Duplicate found for '{sequence_record.id}'.")
            meta_data[sequence_record.id] = sequence_metadata

            sequences_written = write_sequences(
                sequence_record,
                handle
            )

    df = pd.DataFrame(meta_data.values())
    df.to_csv(
        args.output_metadata,
        index=False,
        sep='\t' if args.output_metadata.endswith('tsv') else ','
    )
