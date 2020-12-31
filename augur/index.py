"""
Count occurrence of bases in a set of sequences
"""

from Bio import SeqIO
import Bio.Seq
import Bio.SeqRecord
import csv

def register_arguments(parser):
    parser.add_argument('--sequences', '-s', required=True, help="sequences in fasta or VCF format")
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
    pass


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
    pass


def run(args):
    '''
    counts the number of a,c,t,g,n and other IUPAC nucleotides in a set of sequences
    '''
    #read in files
    try:
        seqs = SeqIO.parse(args.sequences, 'fasta')
    except ValueError as error:
        print("ERROR: Problem reading in {}:".format(args.sequences))
        print(error)
        return 1

    good_chars = {'a', 'c', 'g', 't', 'n', '-', 'r', 'y', 's', 'w', 'k', 'm', 'd', 'h', 'b', 'v','?'}
    tot_length = 0
    num_of_seqs = 0
    #characters to count as list of strings
    other_IUPAC = ['r', 'y', 's', 'w', 'k', 'm', 'd', 'h', 'b', 'v']
    values = [['a'],['c'],['g'],['t'],['n'],other_IUPAC]
    #labels for header in tsv output file
    labels = ['A','C','G','T','N','other_IUPAC']

    #print a warning if labels do not match characters to count 
    if len(labels)-len(values)!=0:
        print('Warning: labels do not  match length of values to count')

    #count occurences and write to output file
    with open(args.output, 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter = '\t')

        #write header
        header = ['sequence_name', 'length']+labels
        tsv_writer.writerow(header)

        for record in seqs:
            #check if sequence contains non nucleotide characters
            if len(set(str(record.seq.lower())).difference(good_chars))==0:

                #change sequences to lower case and remove gaps
                seq = record.seq.ungap(gap='-').ungap(gap='?').lower()
                l = len(seq)
                tot_length += l 
                num_of_seqs += 1    

                #count occurrence for each value
                counts = []
                for v in values:
                    counts.append(sum(map(lambda x: seq.count(x), v)))
                #write row for each sequence
                row = [record.id, l]+counts
                tsv_writer.writerow(row)

            else:
                #skip sequences with non nucleotide characters
                print("sequence contains non nucleotide characters, skipping...")

    if num_of_seqs==0:
        print("No valid sequences found.")
    else:
        print("Analysed %i sequences with an average length of %i nucleotides."%(num_of_seqs,int(tot_length/num_of_seqs)))
        
    seqs.close()