"""
Parse delimited fields from FASTA sequence names into a TSV and FASTA file.
"""

from Bio import SeqIO
import pandas as pd

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
    seqs = SeqIO.parse(args.sequences, 'fasta')

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
    with open(args.output_sequences, 'w', encoding='utf-8') as output:
        for seq in seqs:
            fields = map(str.strip, seq.description.split(args.separator))
            tmp_meta = dict(zip(args.fields, fields))

            tmp_name = tmp_meta[strain_key].translate(forbidden_chactacters)
            seq.name = seq.id = tmp_name
            seq.description = ''

            if args.prettify_fields:
                for field in tmp_meta.keys() & args.prettify_fields:
                    if isinstance(tmp_meta[field], str):
                        tmp_meta[field] = prettify(tmp_meta[field], camelCase=(not field.startswith('author')),
                                                    etal='lower' if field.startswith('author') else None)

            # parse dates and convert to a canonical format
            if args.fix_dates and 'date' in tmp_meta:
                tmp_meta['date'] = fix_dates(tmp_meta['date'],
                                            dayfirst=args.fix_dates=='dayfirst')

            if 'strain' in tmp_meta:
                del tmp_meta['strain']
            meta_data[seq.id] = tmp_meta

            SeqIO.write(seq, output, 'fasta')

    df = pd.DataFrame.from_dict(meta_data, orient='index')
    df.to_csv(args.output_metadata, index_label='strain',
              sep='\t' if args.output_metadata.endswith('tsv') else ',')
