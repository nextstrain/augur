from Bio import SeqIO
from pprint import pprint

with open('../../data/dengue_titers.tsv', 'r') as f:
	titerstrains = set([ line.split()[0] for line in f ])
with open('../../data/dengue_titers.tsv', 'r') as f:
	serastrains = set([ line.split()[1] for line in f ])

autologous = titerstrains.intersection(serastrains)
print len(autologous)

strains_with_titers = [s for s in SeqIO.parse(open('../../data/dengue.fasta', 'r'), 'fasta') if s.description.split('|')[0] in autologous ]
SeqIO.write(strains_with_titers, '../../data/control.fasta', 'fasta')

print 'Found %d strains with autologous titers and sequence data.'%len(strains_with_titers)

