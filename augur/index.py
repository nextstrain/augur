"""
Count occurences of bases in a set of sequences
"""

from Bio import SeqIO
import csv

def register_arguments(parser):
    parser.add_argument('--sequences', '-s', required=True, help="sequences in fasta or VCF format")
    parser.add_argument('--output', '-o', help="output file", required=True)


def run(args):
    '''
    counts the number of a,c,t,g,n and non_actgn bases in a set of sequences
    '''
    #read in files
    try:
        seqs = SeqIO.parse(args.sequences, 'fasta')
    except ValueError as error:
        print("ERROR: Problem reading in {}:".format(args.sequences))
        print(error)
        return 1

    #count occurences and write to output file
    with open(args.output, 'wt') as out_file:
    	tsv_writer = csv.writer(out_file, delimiter = '\t')
    	tsv_writer.writerow(['sequence_name', 'length', '#A','#C', '#G', '#T' ,'#N', '#non-ACTGN'])

    	for record in seqs:
    		#change sequences to lower case and remove gaps
    	   	seq = record.seq.ungap(gap='-').lower()
    	   	l = len(seq)

    	   	#count bases
    	   	a = seq.count('a')
    	   	c = seq.count('c')
    	   	g = seq.count('g')
    	   	t = seq.count('t')
    	   	n = seq.count('n')
    	   	non_acgtn = l-sum([a,c,g,t,n])

    	   	#write row for each sequence
    	   	tsv_writer.writerow([record.id,l,a,c,g,t,n,non_acgtn])

    seqs.close()
    out_file.close()