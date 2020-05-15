"""
Align multiple sequences from FASTA.
"""

import os
from shutil import copyfile
import numpy as np
from Bio import AlignIO, SeqIO, Seq, Align
from .utils import run_shell_command, nthreads_value, shquote

class AlignmentError(Exception):
    # TODO: this exception should potentially be renamed and made augur-wide
    # thus allowing any module to raise it and have the message printed & augur
    # exit with code 1
    pass

def register_arguments(parser):
    parser.add_argument('--sequences', '-s', required=True, nargs="+", metavar="FASTA", help="sequences to align")
    parser.add_argument('--output', '-o', default="alignment.fasta", help="output file (default: %(default)s)")
    parser.add_argument('--nthreads', type=nthreads_value, default=1,
                                help="number of threads to use; specifying the value 'auto' will cause the number of available CPU cores on your system, if determinable, to be used")
    parser.add_argument('--method', default='mafft', choices=["mafft"], help="alignment program to use")
    parser.add_argument('--reference-name', metavar="NAME", type=str, help="strip insertions relative to reference sequence; use if the reference is already in the input sequences")
    parser.add_argument('--reference-sequence', metavar="PATH", type=str, help="Add this reference sequence to the dataset & strip insertions relative to this. Use if the reference is NOT already in the input sequences")
    parser.add_argument('--remove-reference', action="store_true", default=False, help="remove reference sequence from the alignment")
    parser.add_argument('--fill-gaps', action="store_true", default=False, help="If gaps represent missing data rather than true indels, replace by N after aligning.")
    parser.add_argument('--existing-alignment', metavar="FASTA", default=False, help="An existing alignment to which the sequences will be added. The ouput alignment will be the same length as this existing alignment.")
    parser.add_argument('--debug', action="store_true", default=False, help="Produce extra files (e.g. pre- and post-aligner files) which can help with debugging poor alignments.")

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
        seqs = read_sequences(*args.sequences)
        existing_aln = read_alignment(args.existing_alignment) if args.existing_alignment else None

        # if we have been given a reference (strain) name, make sure it is present
        ref_name = args.reference_name
        if args.reference_name:
            ensure_reference_strain_present(ref_name, existing_aln, seqs)

        # If given an existing alignment, then add the reference sequence to this if desired (and if it is the same length)
        if existing_aln and args.reference_sequence:
            existing_aln_fname = args.existing_alignment + ".ref.fasta"
            ref_seq = read_reference(args.reference_sequence)
            if len(ref_seq) != existing_aln.get_alignment_length():
                raise AlignmentError("ERROR: Provided existing alignment ({}bp) is not the same length as the reference sequence ({}bp)".format(existing_aln.get_alignment_length(), len(ref_seq)))
            existing_aln.append(ref_seq)
            write_seqs(existing_aln, existing_aln_fname)
            temp_files_to_remove.append(existing_aln_fname)
            ref_name = ref_seq.id
        else:
            existing_aln_fname = args.existing_alignment # may be False

        ## Create a single file of sequences for alignment (or to be added to the alignment).
        ## Add in the reference file to the sequences _if_ we don't have an existing alignment
        if args.reference_sequence and not existing_aln:
            seqs_to_align_fname = args.output+".to_align.fasta"
            ref_seq = read_reference(args.reference_sequence)
            # reference sequence needs to be the first one for auto direction adjustment (auto reverse-complement)
            write_seqs([ref_seq] + list(seqs.values()), seqs_to_align_fname)
            ref_name = ref_seq.id
        elif existing_aln:
            seqs_to_align_fname = args.output+".new_seqs_to_align.fasta"
            seqs = prune_seqs_matching_alignment(seqs, existing_aln)
            write_seqs(list(seqs.values()), seqs_to_align_fname)
        else:
            seqs_to_align_fname = args.output+".to_align.fasta"
            write_seqs(list(seqs.values()), seqs_to_align_fname)
        temp_files_to_remove.append(seqs_to_align_fname)

        check_duplicates(existing_aln, ref_name, seqs)

        # before aligning, make a copy of the data that the aligner receives as input (very useful for debugging purposes)
        if args.debug and not existing_aln:
            copyfile(seqs_to_align_fname, args.output+".pre_aligner.fasta")

        # generate alignment command & run
        cmd = generate_alignment_cmd(args.method, args.nthreads, existing_aln_fname, seqs_to_align_fname, args.output, args.output+".log")
        success = run_shell_command(cmd)
        if not success:
            raise AlignmentError("Error during alignment")

        # after aligning, make a copy of the data that the aligner produced (useful for debugging)
        if args.debug:
            copyfile(args.output, args.output+".post_aligner.fasta")

        # reads the new alignment
        seqs = read_alignment(args.output)

        # convert the aligner output to upper case and remove auto reverse-complement prefix
        prettify_alignment(seqs)

        # if we've specified a reference, strip out all the columns not present in the reference
        # this will overwrite the alignment file
        if ref_name:
            seqs = strip_non_reference(seqs, ref_name, keep_reference=not args.remove_reference)
        if args.fill_gaps:
            make_gaps_ambiguous(seqs)

        # write the modified sequences back to the alignment file
        write_seqs(seqs, args.output)


    except AlignmentError as e:
        print(str(e))
        return 1

    # finally, remove any temporary files
    for fname in temp_files_to_remove:
        os.remove(fname)

#####################################################################################################

def read_sequences(*fnames):
    seqs = {}
    try:
        for fname in fnames:
            for record in SeqIO.parse(fname, 'fasta'):
                if record.name in seqs and record.seq != seqs[record.name].seq:
                    raise AlignmentError("Detected duplicate input strains \"%s\" but the sequences are different." % record.name)
                    # if the same sequence then we can proceed (and we only take one)
                seqs[record.name] = record
    except FileNotFoundError:
        raise AlignmentError("\nCannot read sequences -- make sure the file %s exists and contains sequences in fasta format" % fname)
    except ValueError as error:
        raise AlignmentError("\nERROR: Problem reading in {}: {}".format(fname, str(error)))
    return seqs

def check_arguments(args):
    # Simple error checking related to a reference name/sequence
    if args.reference_name and args.reference_sequence:
        raise AlignmentError("ERROR: You cannot provide both --reference-name and --reference-sequence")
    if args.remove_reference and not (args.reference_name or args.reference_sequence):
        raise AlignmentError("ERROR: You've asked to remove the reference but haven't specified one!")

def read_alignment(fname):
    try:
        return AlignIO.read(fname, 'fasta')
    except Exception as error:
        raise AlignmentError("\nERROR: Problem reading in {}: {}".format(fname, str(error)))

def ensure_reference_strain_present(ref_name, existing_alignment, seqs):
    if existing_alignment:
        if ref_name not in {x.name for x in existing_alignment}:
            raise AlignmentError("ERROR: Specified reference name %s (via --reference-name) is not in the supplied alignment."%ref_name)
    else:
        if ref_name not in seqs:
            raise AlignmentError("ERROR: Specified reference name %s (via --reference-name) is not in the sequence sample."%ref_name)


    # align
    # if args.method=='mafft':
    #     shoutput = shquote(output)
    #     shname = shquote(seq_fname)
    #     cmd = "mafft --reorder --anysymbol --thread %d %s 1> %s 2> %s.log"%(args.nthreads, shname, shoutput, shoutput)
def read_reference(ref_fname):
    if not os.path.isfile(ref_fname):
        raise AlignmentError("ERROR: Cannot read reference sequence."
                             "\n\tmake sure the file \"%s\" exists"%ref_fname)
    try:
        ref_seq = SeqIO.read(ref_fname, 'genbank' if ref_fname.split('.')[-1] in ['gb', 'genbank'] else 'fasta')
    except:
        raise AlignmentError("ERROR: Cannot read reference sequence."
                "\n\tmake sure the file %s contains one sequence in genbank or fasta format"%ref_fname)
    return ref_seq

def generate_alignment_cmd(method, nthreads, existing_aln_fname, seqs_to_align_fname, aln_fname, log_fname):
    if method=='mafft':
        if existing_aln_fname:
            cmd = "mafft --add %s --keeplength --reorder --anysymbol --nomemsave --adjustdirection --thread %d %s 1> %s 2> %s"%(shquote(seqs_to_align_fname), nthreads, shquote(existing_aln_fname), shquote(aln_fname), shquote(log_fname))
        else:
            cmd = "mafft --reorder --anysymbol --nomemsave --adjustdirection --thread %d %s 1> %s 2> %s"%(nthreads, shquote(seqs_to_align_fname), shquote(aln_fname), shquote(log_fname))
        print("\nusing mafft to align via:\n\t" + cmd +
              " \n\n\tKatoh et al, Nucleic Acid Research, vol 30, issue 14"
              "\n\thttps://doi.org/10.1093%2Fnar%2Fgkf436\n")
    else:
        raise AlignmentError('ERROR: alignment method %s not implemented'%method)
    return cmd

def strip_non_reference(aln, reference, keep_reference=False):
    '''
    return sequences that have all insertions relative to the reference
    removed. The alignment is read from file and returned as list of sequences.

    Parameters
    ----------
    aln : MultipleSeqAlign
        Biopython Alignment
    reference : str
        name of reference sequence, assumed to be part of the alignment
    keep_reference : bool, optional
        by default, the reference sequence is removed after stripping
        non-reference sequence. To keep the reference, use keep_reference=True

    Returns
    -------
    list
        list of trimmed sequences, effectively a multiple alignment

    Tests
    -----
    >>> [s.name for s in strip_non_reference(read_alignment("tests/data/align/test_aligned_sequences.fasta"), "with_gaps", keep_reference=False)]
    Trimmed gaps in with_gaps from the alignment
    ['no_gaps', 'some_other_seq', '_R_crick_strand']
    >>> [s.name for s in strip_non_reference(read_alignment("tests/data/align/test_aligned_sequences.fasta"), "with_gaps", keep_reference=True)]
    Trimmed gaps in with_gaps from the alignment
    ['with_gaps', 'no_gaps', 'some_other_seq', '_R_crick_strand']
    >>> [s.name for s in strip_non_reference(read_alignment("tests/data/align/test_aligned_sequences.fasta"), "no_gaps", keep_reference=True)]
    No gaps in alignment to trim (with respect to the reference, no_gaps)
    ['with_gaps', 'no_gaps', 'some_other_seq', '_R_crick_strand']
    >>> [s.name for s in strip_non_reference(read_alignment("tests/data/align/test_aligned_sequences.fasta"), "no_gaps", keep_reference=False)]
    No gaps in alignment to trim (with respect to the reference, no_gaps)
    ['with_gaps', 'some_other_seq', '_R_crick_strand']
    >>> [s.name for s in strip_non_reference(read_alignment("tests/data/align/test_aligned_sequences.fasta"), "missing", keep_reference=False)]
    Traceback (most recent call last):
      ...
    augur.align.AlignmentError: ERROR: reference missing not found in alignment
    '''
    seqs = {s.name:s for s in aln}
    if reference in seqs:
        ref_array = np.array(seqs[reference])
        if "-" not in ref_array:
            print("No gaps in alignment to trim (with respect to the reference, %s)"%reference)
            return [seq for seq in aln if (keep_reference or seq.name != reference)]
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

def prettify_alignment(aln):
    '''
    Converts all bases to uppercase and removes auto reverse-complement prefix (_R_).
    This modifies the alignment in place.

    Parameters
    ----------
    aln : MultipleSeqAlign
        Biopython Alignment
    '''
    for seq in aln:
        seq.seq = seq.seq.upper()
        # AlignIO.read repeats the ID in the name and description field
        if seq.id.startswith("_R_"):
            seq.id = seq.id[3:]
            print("Sequence \"{}\" was reverse-complemented by the alignment program.".format(seq.id))
        if seq.name.startswith("_R_"):
            seq.name = seq.name[3:]
        if seq.description.startswith("_R_"):
            seq.description = seq.description[3:]

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
        _seq = str(seq.seq)
        _seq = _seq.replace('-', 'N')
        seq.seq = Seq.Seq(_seq, alphabet=seq.seq.alphabet)
        

def check_duplicates(*values):
    names = set()
    def add(name):
        if name in names:
            raise AlignmentError("Duplicate strains of \"{}\" detected".format(name))
        names.add(name)

    for sample in values:
        if not sample:
            # allows false-like values (e.g. always provide existing_alignment, allowing
            # the default which is `False`)
            continue
        elif type(sample) == dict:
            for s in sample:
                add(s)
        elif type(sample) == Align.MultipleSeqAlignment:
            for s in sample:
                add(s.name)
        elif type(sample) == str:
            add(sample)
        else:
            raise TypeError()

def write_seqs(seqs, fname):
    """A wrapper around SeqIO.write with error handling"""
    try:
        SeqIO.write(seqs, fname, 'fasta')
    except FileNotFoundError:
        raise AlignmentError('ERROR: Couldn\'t write "{}" -- perhaps the directory doesn\'t exist?'.format(fname))


def prune_seqs_matching_alignment(seqs, aln):
    """
    Return a set of seqs excluding those set via `exclude` & print a warning
    message for each sequence which is exluded.
    """
    ret = {}
    exclude_names = {s.name for s in aln}
    for name, seq in seqs.items():
        if name in exclude_names:
            print("Excluding {} as it is already present in the alignment".format(name))
        else:
            ret[name] = seq
    return ret
