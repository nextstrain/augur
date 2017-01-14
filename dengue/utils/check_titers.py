from Bio import SeqIO
from pprint import pprint

records = open('dengue_titers.tsv', 'r').readlines()
titerstrains = set(sorted([ l.split()[0] for l in records ]))
serastrains = set(sorted([ l.split()[1] for l in records ]))

autologous = titerstrains.intersection(serastrains)
no_autologous = set(titerstrains).symmetric_difference(serastrains)
all_titers_sera = titerstrains.union(serastrains)

fasta = set(sorted([ s.description.split('|')[0] for s in SeqIO.parse(open('dengue.fasta', 'r'), 'fasta')]))

virus_sera_fasta = autologous.intersection(fasta)
missing_in_fasta = autologous.difference(fasta)

print 'Found %d records in fasta, sera, AND viruses'%len(virus_sera_fasta)

print '\n\nFound %d records in viruses OR sera, but no autologous titers.'%len(no_autologous)
pprint(no_autologous)

print '\n\nFound %d records with autologous titers, but no sequence'%len(missing_in_fasta)
pprint(missing_in_fasta)

if len(no_autologous) > 0 or len(missing_in_fasta) > 0:
	print '\n\n Titers and sera'
	pprint(sorted(list(all_titers_sera)))

strains = set(sorted([ l.split()[0] for l in open('dengue_strains.tsv', 'r').readlines() if l.startswith('DENV3') ]))
missing_strains = virus_sera_fasta.difference(strains)

print 'Found %d strains in virus, sera, and fasta, but not in strain count file'%len(missing_strains)
pprint(missing_strains)
