"""
Align multiple sequences from FASTA or VCF.
"""

import os,sys,argparse
from shutil import copyfile
import numpy as np
from Bio import AlignIO, SeqIO, Seq
from .utils import run_shell_command, nthreads_value

def make_gaps_ambiguous(aln):
    '''
    replace all gaps by 'N' in all sequences in the alignment. TreeTime will treat them
    as fully ambiguous and replace then with the most likely state. This modifies the
    alignment in place.

    Parameters
    ----------
    aln : MultipleSeqAlign
        Biopython Alignment
    '''
    for seq in aln:
        seq_array = np.array(seq)
        gaps = seq_array=='-'
        seq_array[gaps]='N'
        seq.seq = Seq.Seq("".join(seq_array))


def register_arguments(parser):
    parser.add_argument('--sequences', '-s', required=True, help="sequences in fasta or VCF format")
    parser.add_argument('--output', '-o', help="output file")
    parser.add_argument('--nthreads', type=nthreads_value, default=1,
                                help="number of threads to use; specifying the value 'auto' will cause the number of available CPU cores on your system, if determinable, to be used")
    parser.add_argument('--method', default='mafft', choices=["mafft"],
                                help="alignment program to use")
    parser.add_argument('--reference-name', type=str, help="strip insertions relative to reference sequence; use if the reference is already in the input sequences")
    parser.add_argument('--reference-sequence', type=str, help="strip insertions relative to reference sequence; use if the reference is NOT already in the input sequences")
    parser.add_argument('--remove-reference', action="store_true", default=False, help="remove reference sequence from the alignment")
    parser.add_argument('--fill-gaps', action="store_true", default=False, help="if gaps represent missing data rather than true indels, replace by N after aligning")


def run(args):
    '''
    Parameters
    ----------
    args : namespace
        arguments passed in via the command-line from augur

    Returns
    -------
    int
        returns 0 for success, 1 for general error
    '''
    seq_fname = args.sequences
    ref_name = args.reference_name
    ref_fname = args.reference_sequence
    if args.output:
        output = args.output
    else:
        output = "aligment.fasta"

    try:
        seqs = {s.id:s for s in SeqIO.parse(seq_fname, 'fasta')}
    except:
        print("Cannot read sequences -- make sure the file %s exists and contains sequences in fasta format"%seq_fname)
        return 1

    if ref_name and (ref_name not in seqs):
        print("Specified reference name %s is not in the sequence sample. Will not trim."%ref_name)
        ref_name = None

    # potentially add the reference seuqnce to the sequences
    if ref_fname and (not ref_name):
        if os.path.isfile(ref_fname):
            try:
                ref_seq = SeqIO.read(ref_fname, 'genbank' if ref_fname.split('.')[-1] in ['gb', 'genbank'] else 'fasta')
            except:
                print("WARNING: Cannot read reference sequence."
                      "\n\tmake sure the file %s contains one sequence in genbank or fasta format"%ref_fname)
            else:
                ref_name = ref_seq.id
                seq_fname+='.ref.fasta'
                SeqIO.write(list(seqs.values())+[ref_seq], seq_fname, 'fasta')
        else:
            print("WARNING: Cannot read reference sequence."
                  "\n\tmake sure the file %s does not exist"%ref_fname)

    # before aligning, make a copy of the data that the aligner receives as input (very useful for debugging purposes)
    copyfile(seq_fname, output+".pre_aligner.fasta")

    # align
    if args.method=='mafft':
        cmd = "mafft --reorder --anysymbol --thread %d %s 1> %s 2> %s.log"%(args.nthreads, seq_fname, output, output)
        print("\nusing mafft to align via:\n\t" + cmd +
              " \n\n\tKatoh et al, Nucleic Acid Research, vol 30, issue 14"
              "\n\thttps://doi.org/10.1093%2Fnar%2Fgkf436\n")
    else:
        print('ERROR: alignment method not implemented')
        return 1

    success = run_shell_command(cmd)
    if not success: # return error if aligner errored
        return 1

    # after aligning, make a copy of the data that the aligner produced (useful for debugging)
    copyfile(output, output+".post_aligner.fasta")

    # convert the aligner output to upper case (replacing the file in place)
    aln = AlignIO.read(output, 'fasta')
    for seq in aln:
        seq.seq = seq.seq.upper()
    AlignIO.write(aln, output, 'fasta')

    # if there's a valid reference in the alignment, strip out all the columns not present in the reference
    # this will overwrite the alignment file
    if ref_name:
        seqs = strip_non_reference(output, ref_name, keep_reference=not args.remove_reference)
        if not seqs:
            return # error already printed from strip_non_reference
        if args.fill_gaps:
            make_gaps_ambiguous(seqs)
        SeqIO.write(seqs, output, 'fasta')


def strip_non_reference(alignment_fname, reference, keep_reference=False):
    '''
    return sequences that have all insertions relative to the reference
    removed. The alignment is read from file and returned as list of sequences.

    Parameters
    ----------
    alignment_fname : str
        alignment file name, file needs to be fasta format
    reference : str
        name of reference sequence, assumed to be part of the alignment
    keep_reference : bool, optional
        by default, the reference sequence is removed after stripping
        non-reference sequence. To keep the reference, use keep_reference=True

    Returns
    -------
    list
        list of trimmed sequences, effectively a multiple alignment
    '''
    aln = AlignIO.read(alignment_fname, 'fasta')
    seqs = {s.name:s for s in aln}
    if reference in seqs:
        ref_array = np.array(seqs[reference])
        ungapped = ref_array!='-'
        ref_aln_array = np.array(aln)[:,ungapped]
    else:
        print("WARNING: reference", reference, "not found in alignment")
        return

    out_seqs = []
    for seq, seq_array in zip(aln, ref_aln_array):
        seq.seq = Seq.Seq(''.join(seq_array))
        if keep_reference or seq.name!=reference:
            out_seqs.append(seq)

    print("Trimmed gaps in", reference, "from the alignment")
    return out_seqs
