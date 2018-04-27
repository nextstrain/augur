import os,sys,argparse
import numpy as np
from Bio import AlignIO, SeqIO, Seq


def run(args):
    parser = argparse.ArgumentParser("Align sequences")
    parser.add_argument('-s', required=True, help="sequences in fasta format")
    parser.add_argument('-o', required=True, help="output file")
    parser.add_argument('--nthreads', type=int, default=2,
                        help="number of threads used by mafft")
    parser.add_argument('--aligner', default='mafft',
                        help="alignment program to use")
    parser.add_argument('--reference', type=str, help="strip insertions relative to reference sequence")
    parser.add_argument('--remove_reference', action="store_true", help="keep reference sequence in alignment")
    args = parser.parse_args(args)

    if args.aligner=='mafft':
        os.system("mafft --anysymbol --thread %d %s 1> %s 2>mafft_stderr"%(args.nthreads, args.s, args.o))
    else:
        print('not implemented')

    from Bio import AlignIO
    aln = AlignIO.read(args.o, 'fasta')
    aln_dict = {}
    for seq in aln:
        seq.seq = seq.seq.upper()
        aln_dict[seq.name]=seq
    AlignIO.write(aln, args.o, 'fasta')


    if args.reference:
        seqs = strip_non_reference(args.o, args.reference, keep_reference=not args.remove_reference)
        SeqIO.write(seqs, args.o, 'fasta')



def strip_non_reference(alignment_fname, reference, keep_reference=False):
    '''
    return sequences that have all insertions relative to the reference
    removed. The alignment is read from file and returned as list of sequences.
    '''
    aln = AlignIO.read(alignment_fname, 'fasta')
    seqs = {s.name:s for s in aln}
    if reference in seqs:
        ref_array = np.array(seqs[reference])
        ungapped = ref_array!='-'
        ref_aln_array = np.array(aln)[:,ungapped]
    else:
        print("reference", reference, "not found in alignment")
        return

    out_seqs = []
    for seq, seq_array in zip(aln, ref_aln_array):
        seq.seq = Seq.Seq(''.join(seq_array))
        if keep_reference or seq.name!=reference:
            out_seqs.append(seq)

    return out_seqs


if __name__ == '__main__':
    run(sys.argv)