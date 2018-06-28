from Bio import SeqIO


"""
Compile a list of sequence names for all egg-passaged sequences and any cell-passed and/or unpassaged matches for these strains
"""

egg_seqs = []
force_include = []
cell_matches = 0
unpassaged_matches = 0
both_matches = 0

#Find all egg-passaged sequences
for seq_record in SeqIO.parse('../../../fauna/data/h3n2_ha.fasta', "fasta"):
    virus_name = str.split(seq_record.id, '|')[0]

    if 'egg' in virus_name:
        egg_seqs.append(virus_name)

#Find all pairs for egg-passaged sequences and add to list of viruses to include during augur sampling
for egg_seq in egg_seqs:
    seq_name = str.split(egg_seq,'-egg')[0]
    cell=0
    unpassaged=0

    for seq_record in SeqIO.parse('../../../fauna/data/h3n2_ha.fasta', "fasta"):
        if seq_name in seq_record.id:
            force_include.append(str.split(seq_record.id, '|')[0])

            if '-egg' not in seq_record.id:
                if '-cell' in seq_record.id:
                    cell+=1
                else:
                    unpassaged+=1

    cell_matches+=cell
    unpassaged_matches+=unpassaged
    if cell!=0 and unpassaged!=0:
        both_matches+=1

#print summary stats
print('Egg-passaged sequences in data file: ' + str(len(egg_seqs)))
print('Egg seqs with matched cell-passaged strain: ' + str(cell_matches))
print('Egg seqs with matched unpassaged strain: ' + str(unpassaged_matches))
print('Egg seqs with matched cell-passaged and unpassaged strains: ' + str(both_matches))

#write force_include list to file
with open('include_seqs.txt', 'w') as filehandle:
    filehandle.writelines('%s\n' %strain_name for strain_name in force_include)
