"""
Align multiple sequences from FASTA.
"""

import os,sys,argparse
from shutil import copyfile
import numpy as np
from Bio import AlignIO, SeqIO, Seq
from .utils import run_shell_command, nthreads_value

class AlignmentError(Exception):
    # TODO: this exception should potentially be renamed and made augur-wide
    # thus allowing any module to raise it and have the message printed & augur
    # exit with code 1
    pass

def register_arguments(parser):
    parser.add_argument('--sequences', '-s', required=True, metavar="FASTA", help="sequences to align")
    parser.add_argument('--output', '-o', default="alignment.fasta", help="output file (default: %(default)s)")
    parser.add_argument('--nthreads', type=nthreads_value, default=1,
                                help="number of threads to use; specifying the value 'auto' will cause the number of available CPU cores on your system, if determinable, to be used")
    parser.add_argument('--method', default='mafft', choices=["mafft"], help="alignment program to use")
    parser.add_argument('--reference-name', metavar="NAME", type=str, help="strip insertions relative to reference sequence; use if the reference is already in the input sequences")
    parser.add_argument('--reference-sequence', metavar="PATH", type=str, help="Add this reference sequence to the dataset & strip insertions relative to this. Use if the reference is NOT already in the input sequences")
    parser.add_argument('--remove-reference', action="store_true", default=False, help="remove reference sequence from the alignment")
    parser.add_argument('--fill-gaps', action="store_true", default=False, help="if gaps represent missing data rather than true indels, replace by N after aligning. A reference must be specified.")


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
    temp_files_to_remove = []

    try:
        check_arguments(args)
        seqs = read_sequences(args.sequences)

        if args.reference_name and (args.reference_name not in seqs):
            raise AlignmentError("ERROR: Specified reference name %s (via --reference-name) is not in the sequence sample."%args.reference_name)
        else:
            ref_name = args.reference_name


        # potentially add the reference sequence to the sequences
        if args.reference_sequence:
            seq_fname = args.sequences+".ref.fasta"
            ref_name = add_reference_seq(args.reference_sequence, seqs, seq_fname)
            temp_files_to_remove.append(seq_fname) # remove this new file upon success
        else:
            seq_fname = args.sequences

        # before aligning, make a copy of the data that the aligner receives as input (very useful for debugging purposes)
        copyfile(seq_fname, args.output+".pre_aligner.fasta")

        # generate alignment command & run
        cmd = generate_alignment_cmd(args.method, args.nthreads, seq_fname, args.output, args.output+".log")
        success = run_shell_command(cmd)
        if not success:
            raise AlignmentError("Error during alignment")

        # after aligning, make a copy of the data that the aligner produced (useful for debugging)
        copyfile(args.output, args.output+".post_aligner.fasta")

        # convert the aligner output to upper case (replacing the file in place)
        write_uppercase_alignment_in_place(args.output)

        # if we've specified a reference, strip out all the columns not present in the reference
        # this will overwrite the alignment file
        if ref_name:
            seqs = strip_non_reference(args.output, ref_name, keep_reference=not args.remove_reference)
            if args.fill_gaps:
                make_gaps_ambiguous(seqs)
            SeqIO.write(seqs, args.output, 'fasta')


    except AlignmentError as e:
        print(str(e))
        return 1

    # finally, remove any temporary files
    if temp_files_to_remove:
        for fname in temp_files_to_remove:
            os.remove(fname)




def read_sequences(fname):
    try:
        return SeqIO.to_dict(SeqIO.parse(fname, 'fasta'))
    except FileNotFoundError:
        raise AlignmentError("\nCannot read sequences -- make sure the file %s exists and contains sequences in fasta format"%fname)
    except ValueError as error:
        raise AlignmentError("\nERROR: Problem reading in {}: {}".format(fname, str(error)))

def check_arguments(args):
    # Simple error checking related to a reference name/sequence
    if args.reference_name and args.reference_sequence:
        raise AlignmentError("ERROR: You cannot provide both --reference-name and --reference-sequence")
    if args.remove_reference and not (args.reference_name or args.reference_sequence):
        raise AlignmentError("ERROR: You've asked to remove the reference but haven't specified one!")
    if args.fill_gaps and not (args.reference_name or args.reference_sequence):
        raise AlignmentError("ERROR: In order to fill gaps (--fill-gaps) you must specify a reference")

def add_reference_seq(ref_fname, seqs, combined_fname):
    # add the reference sequence found in ref_fname to the sequences in seqs
    # and write out to the filename combined_fname
    if not os.path.isfile(ref_fname):
        raise AlignmentError("ERROR: Cannot read reference sequence."
                             "\n\tmake sure the file \"%s\" exists"%ref_fname)
    try:
        ref_seq = SeqIO.read(ref_fname, 'genbank' if ref_fname.split('.')[-1] in ['gb', 'genbank'] else 'fasta')
    except:
        raise AlignmentError("ERROR: Cannot read reference sequence."
                "\n\tmake sure the file %s contains one sequence in genbank or fasta format"%ref_fname)
    SeqIO.write(list(seqs.values())+[ref_seq], combined_fname, 'fasta')
    return ref_seq.id


def generate_alignment_cmd(method, nthreads, seqs_to_align, aln_fname, log_fname):
    if method=='mafft':
        cmd = "mafft --reorder --anysymbol --thread %d %s 1> %s 2> %s"%(nthreads, seqs_to_align, aln_fname, log_fname)
        print("\nusing mafft to align via:\n\t" + cmd +
              " \n\n\tKatoh et al, Nucleic Acid Research, vol 30, issue 14"
              "\n\thttps://doi.org/10.1093%2Fnar%2Fgkf436\n")
    else:
        raise AlignmentError('ERROR: alignment method %s not implemented'%method)
    return cmd


def write_uppercase_alignment_in_place(fname):
    aln = AlignIO.read(fname, 'fasta')
    for seq in aln:
        seq.seq = seq.seq.upper()
    AlignIO.write(aln, fname, 'fasta')



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
        raise AlignmentError("ERROR: reference %s not found in alignment"%reference)

    out_seqs = []
    for seq, seq_array in zip(aln, ref_aln_array):
        seq.seq = Seq.Seq(''.join(seq_array))
        if keep_reference or seq.name!=reference:
            out_seqs.append(seq)

    print("Trimmed gaps in", reference, "from the alignment")
    return out_seqs


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