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


def register_arguments(parser):
    parser.add_argument('--sequences', '-s', required=True, help="sequences in fasta or VCF format")
    parser.add_argument('--output-sequences', help="output sequences file")
    parser.add_argument('--output-metadata', help="output metadata file")
    parser.add_argument('--fields', nargs='+', help="fields in fasta header")
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
        fields = seq.description.split(args.separator)
        tmp_name = fields[strain_index]
        for x, y in forbidden_characters:
            tmp_name = tmp_name.replace(x, y)

        seq.name = seq.id = tmp_name
        seq.description = ''
        meta_data[seq.id] = {k:v for k,v in zip(args.fields, fields) }
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
