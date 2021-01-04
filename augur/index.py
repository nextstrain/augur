"""
Count occurrence of bases in a set of sequences
"""

from Bio import SeqIO
import Bio.Seq
import Bio.SeqRecord
import sys
import csv

def register_arguments(parser):
    parser.add_argument('--sequences', '-s', required=True, help="sequences in fasta format")
    parser.add_argument('--output', '-o', help="output file", required=True)


def index_sequence(sequence, values):
    """Count the number of nucleotides for a given sequence record.

    Parameters
    ----------
    sequence : Bio.SeqRecord.SeqRecord
        sequence record to index.

    values : list of lists of str
        values to count

    Returns
    -------
    list :
        summary of the given sequence's strain name, length, nucleotide counts
        for the given values, and a final column with the number of characters
        that didn't match any of those in the given values.

    >>> other_IUPAC = ['r', 'y', 's', 'w', 'k', 'm', 'd', 'h', 'b', 'v']
    >>> values = [['a'],['c'],['g'],['t'],['n'], other_IUPAC, ['-'], ['?']]
    >>> sequence_a = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq("ACTGN-?XWN"), id="seq_A")
    >>> index_sequence(sequence_a, values)
    ['seq_A', 10, 1, 1, 1, 1, 2, 2, 1, 1, 0]

    >>> sequence_b = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq("ACTGACTG"), id="seq_B")
    >>> index_sequence(sequence_b, values)
    ['seq_B', 8, 2, 2, 2, 2, 0, 0, 0, 0, 0]

    Characters in the given sequence that are not in the given list of values to
    count get counted in the final column.

    >>> sequence_c = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq("ACTG%@!!!NN"), id="seq_C")
    >>> index_sequence(sequence_c, values)
    ['seq_C', 11, 1, 1, 1, 1, 1, 0, 0, 0, 5]

    """
    counts = []
    seq = sequence.seq.lower()
    l = len(seq)   
    
    for v in values:
        counts.append(sum(map(lambda x: seq.count(x), v)))
                
    invalid_nucleotides = l-sum(counts)
    row = [sequence.id, l]+counts+[invalid_nucleotides]    
    
    return row


def index_sequences(sequences_path, sequence_index_path):
    """Count the number of A, C, T, G, N, other IUPAC nucleotides, ambiguous bases
    ("?" and "-"), and other invalid characters in a set of sequences and write
    the composition as a data frame to the given sequence index path.

    Parameters
    ----------
    sequences_path : str or Path-like
        path to a sequence file to index.

    sequence_index_path : str or Path-like
        path to a tab-delimited file containing the composition details for each
        sequence in the given input file.

    """
    #read in files
    try:
        seqs = SeqIO.parse(sequences_path, 'fasta')
    except ValueError as error:
        print("ERROR: Problem reading in {}:".format(sequences_path), file=sys.stderr)
        print(error, file=sys.stderr)
        return 1
    
    other_IUPAC = ['r', 'y', 's', 'w', 'k', 'm', 'd', 'h', 'b', 'v']
    values = [['a'],['c'],['g'],['t'],['n'],other_IUPAC,['-'],['?']] 
    labels = ['A','C','G','T','N','other_IUPAC','-','?']
    
    tot_length = 0
    num_of_seqs = 0
    
    with open(sequence_index_path, 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter = '\t')
        
        #write header i output file
        header = ['strain', 'length']+labels+['invalid_nucleotides']
        tsv_writer.writerow(header)
        
        for record in seqs:
            #index the sequence and write row in output file
            row = index_sequence(record, values)
            tsv_writer.writerow(row)
            
            tot_length += row[1]
            num_of_seqs += 1 
    
    seqs.close()    
    print("Analysed %i sequences with an average length of %i nucleotides."%(num_of_seqs,int(tot_length/num_of_seqs)))


def run(args):
    '''
    runs index_sequences which counts the number of A, C, T, G, N, other IUPAC nucleotides, ambiguous bases
    ("?" and "-"), and other invalid characters in a set of sequences and write
    the composition as a data frame to the given sequence index path.
    '''
    print('Running index_sequences:')
    index_sequences(args.sequences, args.output)