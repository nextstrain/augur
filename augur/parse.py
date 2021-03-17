"""
Parse delimited fields from FASTA sequence names into a TSV and FASTA file.
"""
import pandas as pd

from .io import open_file, read_sequences, write_sequences

forbidden_chactacters = str.maketrans(
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

def fix_dates(d, dayfirst=True):
    '''
    attempt to parse a date string using pandas date parser. If ambiguous,
    the argument 'dayfirst' determines whether month or day is assumed to be
    the first field. Incomplete dates will be padded with XX.
    On failure to parse the date, the function will return the input.
    '''
    try:
        from pandas.core.tools.datetimes import parsing
        results = parsing.parse_time_string(d, dayfirst=dayfirst)
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
    except Exception as e:
        print("WARNING: unable to parse %s as date"%d, e)
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


def parse_sequence(sequence, fields, strain_key="strain", separator="|", prettify_fields=None, fix_dates=None):
    """Parse a single sequence record into a sequence record and associated metadata.

    Parameters
    ----------
    sequence : Bio.SeqRecord.SeqRecord
        a BioPython sequence record to parse with metadata stored in its description field.

    fields : list or tuple
        a list of names for fields expected in the given record's description.

    strain_key : str
        name of the field to use as the given sequence's unique id

    separator : str
        delimiter to split record description by.

    prettify_fields : list or tuple
        a list of field names for which the values in those fields should be prettified.

    fix_dates : str
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

    tmp_name = metadata[strain_key].translate(forbidden_chactacters)
    sequence.name = sequence.id = tmp_name
    sequence.description = ''

    if prettify_fields:
        for field in metadata.keys() & prettify_fields:
            if isinstance(metadata[field], str):
                metadata[field] = prettify(metadata[field], camelCase=(not field.startswith('author')),
                                            etal='lower' if field.startswith('author') else None)

    # parse dates and convert to a canonical format
    if fix_dates and 'date' in metadata:
        metadata['date'] = fix_dates(
            metadata['date'],
            dayfirst=fix_dates=='dayfirst'
        )

    metadata["strain"] = sequence.id

    return sequence, metadata


def register_arguments(parser):
    parser.add_argument('--sequences', '-s', required=True, help="sequences in fasta or VCF format")
    parser.add_argument('--output-sequences', help="output sequences file")
    parser.add_argument('--output-metadata', help="output metadata file")
    parser.add_argument('--fields', nargs='+', help="fields in fasta header")
    parser.add_argument('--prettify-fields', nargs='+', help="apply string prettifying operations (underscores to spaces, capitalization, etc) to specified metadata fields")
    parser.add_argument('--separator', default='|', help="separator of fasta header")
    parser.add_argument('--fix-dates', choices=['dayfirst', 'monthfirst'],
                                help="attempt to parse non-standard dates and output them in standard YYYY-MM-DD format")


def run(args):
    '''
    parse a fasta file and turn information in the header into
    a tsv or csv file.
    '''
    sequences = read_sequences(args.sequences)

    # if strain or name are found in specified fields, use this
    # field to index the dictionary and the data frame
    meta_data = {}

    if 'name' in args.fields:
        strain_key = 'name'
    elif 'strain' in args.fields:
        strain_key = 'strain'
    else:
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
