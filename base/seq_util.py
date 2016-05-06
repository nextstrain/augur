nuc_alpha = 'ACGT-N'
aa_alpha = 'ACDEFGHIKLMNPQRSTVWY*-X'

def pad_nucleotide_sequences(aln_aa, seq_nuc):
    '''
    introduce gaps of 3 (---) into nucleotide sequences corresponding to aligned DNA sequences.

    Parameters:
    - aln_aa: amino acid alignment
    - seq_nuc: unaligned nucleotide sequences.

    Returns:
    - aligned nucleotide sequences with all gaps length 3
    '''
    from Bio.Align import MultipleSeqAlignment
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    aln_nuc = MultipleSeqAlignment([])
    for aa_seq  in aln_aa:
        try:
            tmp_nuc_seq = str(seq_nuc[aa_seq.id].seq)
        except KeyError as e:
            print aa_seq.id
            print 'Key not found, continue with next sequence'
            continue

        tmpseq = ''
        nuc_pos = 0
        for aa in aa_seq:
            if aa=='-':
                tmpseq+='---'
            else:
                tmpseq+=tmp_nuc_seq[nuc_pos:(nuc_pos+3)]
                nuc_pos+=3

        aln_nuc.append(SeqRecord(seq=Seq(tmpseq),id=aa_seq.id))

    return aln_nuc
