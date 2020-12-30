"""
Count occurrence of bases in a set of sequences
"""

from Bio import SeqIO
import csv

def register_arguments(parser):
    parser.add_argument('--sequences', '-s', required=True, help="sequences in fasta or VCF format")
    parser.add_argument('--output', '-o', help="output file", required=True)


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