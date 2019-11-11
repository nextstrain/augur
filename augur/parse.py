"""
Parse delimited fields from FASTA sequence names into a TSV and FASTA file.
"""

from Bio import SeqIO
import pandas as pd

forbidden_characters = [
    (' ',''),
    ('(','_'),
    (')','_'),
    ('[','_'),
    (']','_'),
    (':','_'),
    (',','_'),
    (';','_'),
    ('\\','_')
]

def fix_dates(d, dayfirst=True):
    '''
    attempt to parse a date string using pandas date parser. If ambiguous,
    the argument 'dayfirst' determines whether month or day is assumed to be
    the first field. Incomplete dates will be padded with XX.
    On failure to parse the date, the function will return the input.
    '''
    try:
        dto, _, res = pd.core.tools.datetimes.parse_time_string(d, dayfirst=dayfirst)
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
    if (trim > 0 and len(x) > trim):
        res = x[:trim] + "..."

    if any(c == res for c in ["usvi", "usa", "uk"]):
        res = res.upper()

    words = res.split('_')

    if (camelCase):
        words = [w[0].upper()+w[1:] for w in words if len(w)]

    res = ' '.join(words)

    if removeComma:
        res.replace(',', '')

    if etal=='lower':
        res = res.replace('Et Al', 'et al')
    elif etal=='strip':
        res = res.replace('et al.', '').replace('Et Al.', '').replace('et al', '').replace('Et Al', '');

    return res;


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
    seqs = list(SeqIO.parse(args.sequences, 'fasta'))

    # if strain or name are found in specified fields, use this
    # field to index the dictionary and the data frame
    meta_data = {}
    if 'name' in args.fields:
        strain_index = args.fields.index('name')
    elif 'strain' in args.fields:
        strain_index = args.fields.index('strain')
    else:
        strain_index = 0

    # loop over sequences, parse fasta header of each sequence
    for seq in seqs:
        fields = [x.strip() for x in seq.description.split(args.separator)]
        tmp_name = fields[strain_index]
        for x, y in forbidden_characters:
            tmp_name = tmp_name.replace(x, y)

        seq.name = seq.id = tmp_name
        seq.description = ''
        tmp_meta = {k:v for k,v in zip(args.fields, fields)}

        if args.prettify_fields:
            for field in args.prettify_fields:
                if field in tmp_meta and type(tmp_meta[field])==str:
                    tmp_meta[field] = prettify(tmp_meta[field], camelCase=(not field.startswith('author')),
                                                  etal='lower' if field.startswith('author') else None)

        meta_data[seq.id] = tmp_meta
        meta_data[seq.id].pop('strain')
        # parse dates and convert to a canonical format
        if args.fix_dates and 'date' in args.fields:
            meta_data[seq.id]['date'] = fix_dates(meta_data[seq.id]['date'],
                                        dayfirst=args.fix_dates=='dayfirst')

    # output results to a new fasta alignment and tsv/csv via pandas
    SeqIO.write(seqs, args.output_sequences, 'fasta')
    df = pd.DataFrame.from_dict(meta_data, orient='index')
    df.to_csv(args.output_metadata, index_label='strain',
              sep='\t' if args.output_metadata[-3:]=='tsv' else ',')
