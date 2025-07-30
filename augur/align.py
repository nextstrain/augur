"""
Align multiple nucleotide sequences from FASTA. The "N" character is treated as
missing or ambiguous sites, so aligning amino acid sequences is not supported.
"""

import os
from shlex import quote as shquote
from shutil import copyfile
import numpy as np
from Bio import AlignIO, SeqIO, Seq, Align
from .argparse_ import ExtendOverwriteDefault
from .io.file import open_file
from .io.sequences import read_sequences, read_single_sequence
from .io.shell_command_runner import run_shell_command
from .utils import nthreads_value
from collections import defaultdict

class AlignmentError(Exception):
    # TODO: this exception should potentially be renamed and made augur-wide
    # thus allowing any module to raise it and have the message printed & augur
    # exit with code 1
    pass

def register_arguments(parser):
    """
    Add arguments to parser.
    Kept as a separate function than `register_parser` to continue to support
    unit tests that use this function to create argparser.
    """
    parser.add_argument('--sequences', '-s', required=True, nargs="+", action=ExtendOverwriteDefault, metavar="FASTA", help="sequences to align")
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


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("align", help=__doc__)
    register_arguments(parser)
    return parser


def prepare(sequences, existing_aln_fname, output, ref_name, ref_seq_fname):
    """Prepare the sequences, existing alignment, and reference sequence for alignment.

    This function:
        1. Combines all given input sequences into a single file
        2. Checks to make sure the input sequences don't overlap with the existing alignment, if one exists.
        3. If given a reference name, check that sequence exists in either the existing alignment, if given, or the input sequences.
        4. If given a reference sequence, either add it to the existing alignment or prepend it to the input seqeunces.
        5. Write the input sequences to a single file, and write the alignment back out if we added the reference sequence to it.

    Parameters
    ----------
    sequences : list of str
        List of paths to FASTA-formatted sequences to align.
    existing_aln_fname : str
        Path of an existing alignment to use, or None
    output: str
        Path the aligned sequences will be written out to.
    ref_name: str
        The name of the reference sequence, if provided
    ref_seq_fname: str
        The path to the reference sequence file. If this is provided, it overrides ref_name.

    Returns
    -------
    tuple of str
        The existing alignment filename, the new sequences filename, and the name of the reference sequence.
    """
    seqs = read_and_validate_sequences(*sequences)
    seqs_to_align_fname = output + ".to_align.fasta"

    if existing_aln_fname:
        existing_aln = read_alignment(existing_aln_fname)
        seqs = prune_seqs_matching_alignment(seqs, existing_aln)
    else:
        existing_aln = None

    if ref_seq_fname:
        ref_seq = read_reference(ref_seq_fname)
        ref_name = ref_seq.id
        if existing_aln:
            if len(ref_seq) != existing_aln.get_alignment_length():
                raise AlignmentError("ERROR: Provided existing alignment ({}bp) is not the same length as the reference sequence ({}bp)".format(existing_aln.get_alignment_length(), len(ref_seq)))
            existing_aln_fname = existing_aln_fname + ".ref.fasta"
            existing_aln.append(ref_seq)
            write_seqs(existing_aln, existing_aln_fname)
        else:
            # reference sequence needs to be the first one for auto direction
            # adjustment (auto reverse-complement)
            seqs.insert(0, ref_seq)
    elif ref_name:
        ensure_reference_strain_present(ref_name, existing_aln, seqs)

    write_seqs(seqs, seqs_to_align_fname)

    # 90% sure this is only ever going to catch ref_seq was a dupe
    check_duplicates(existing_aln, seqs)
    return existing_aln_fname, seqs_to_align_fname, ref_name

def run(args):
    '''
    Parameters
    ----------
    args : argparse.Namespace
        arguments passed in via the command-line from augur

    Returns
    -------
    int
        returns 0 for success, 1 for general error
    '''
    temp_files_to_remove = []

    try:
        check_arguments(args)
        existing_aln_fname, seqs_to_align_fname, ref_name = prepare(args.sequences, args.existing_alignment, args.output, args.reference_name, args.reference_sequence)
        temp_files_to_remove.append(seqs_to_align_fname)
        if existing_aln_fname != args.existing_alignment:
            temp_files_to_remove.append(existing_aln_fname)
        # -- existing_aln_fname, seqs_to_align_fname, ref_name --

        # before aligning, make a copy of the data that the aligner receives as input (very useful for debugging purposes)
        if args.debug and not existing_aln_fname:
            copyfile(seqs_to_align_fname, args.output+".pre_aligner.fasta")

        # generate alignment command & run
        log = args.output + ".log"
        cmd = generate_alignment_cmd(args.method, args.nthreads, existing_aln_fname, seqs_to_align_fname, args.output, log)
        success = run_shell_command(cmd)
        if not success:
            raise AlignmentError(f"Error during alignment: please see the log file {log!r} for more details")

        # after aligning, make a copy of the data that the aligner produced (useful for debugging)
        if args.debug:
            copyfile(args.output, args.output+".post_aligner.fasta")

        postprocess(args.output, ref_name, not args.remove_reference, args.fill_gaps)


    except AlignmentError as e:
        print(str(e))
        return 1

    # finally, remove any temporary files
    for fname in temp_files_to_remove:
        os.remove(fname)


def postprocess(output_file, ref_name, keep_reference, fill_gaps):
    """Postprocessing of the combined alignment file.

    The modified alignment is written directly to output_file.

    Parameters
    ----------
    output_file: str
        The file the new alignment was written to
    ref_name: str
        If provided, the name of the reference strain used in the alignment
    keep_reference: bool
        If the reference was provided, whether it should be kept in the alignment
    fill_gaps: bool
        Replace all gaps in the alignment with "N" to indicate ambiguous sites.
    """
    # -- ref_name --
    # reads the new alignment
    seqs = read_alignment(output_file)
    # convert the aligner output to upper case and remove auto reverse-complement prefix
    prettify_alignment(seqs)

    # if we've specified a reference, strip out all the columns not present in the reference
    # this will overwrite the alignment file
    if ref_name:
        seqs = strip_non_reference(seqs, ref_name, insertion_csv=output_file+".insertions.csv")
        if not keep_reference:
            seqs = remove_reference_sequence(seqs, ref_name)

    if fill_gaps:
        make_gaps_ambiguous(seqs)

    # write the modified sequences back to the alignment file
    write_seqs(seqs, output_file)



#####################################################################################################

# Note: This is effectively augur.io.read_sequences with extra validation.
def read_and_validate_sequences(*fnames):
    """return list of sequences from all fnames"""
    seqs = {}
    try:
        for fname in fnames:
            for record in read_sequences(fname, format='fasta'):
                if record.name in seqs and record.seq != seqs[record.name].seq:
                    raise AlignmentError("Detected duplicate input strains \"%s\" but the sequences are different." % record.name)
                    # if the same sequence then we can proceed (and we only take one)
                seqs[record.name] = record
    except FileNotFoundError:
        raise AlignmentError("\nCannot read sequences -- make sure the file %s exists and contains sequences in fasta format" % fname)
    except ValueError as error:
        raise AlignmentError("\nERROR: Problem reading in {}: {}".format(fname, str(error)))
    return list(seqs.values())

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
        if ref_name not in {x.name for x in seqs}:
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
        ref_seq = read_single_sequence(ref_fname, format='genbank' if ref_fname.split('.')[-1] in ['gb', 'genbank'] else 'fasta')
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


def remove_reference_sequence(seqs, reference_name):
    return [seq for seq in seqs if seq.name!=reference_name]


def strip_non_reference(aln, reference, insertion_csv=None):
    '''
    return sequences that have all insertions relative to the reference
    removed. The aligment is returned as list of sequences.

    Parameters
    ----------
    aln : Bio.Align.MultipleSeqAlignment
        Biopython Alignment
    reference : str
        name of reference sequence, assumed to be part of the alignment

    Returns
    -------
    list
        list of trimmed sequences, effectively a multiple alignment

    Examples
    --------
    >>> [s.name for s in strip_non_reference(read_alignment("tests/data/align/test_aligned_sequences.fasta"), "with_gaps")]
    Trimmed gaps in with_gaps from the alignment
    ['with_gaps', 'no_gaps', 'some_other_seq', '_R_crick_strand']
    >>> [s.name for s in strip_non_reference(read_alignment("tests/data/align/test_aligned_sequences.fasta"), "no_gaps")]
    No gaps in alignment to trim (with respect to the reference, no_gaps)
    ['with_gaps', 'no_gaps', 'some_other_seq', '_R_crick_strand']
    >>> [s.name for s in strip_non_reference(read_alignment("tests/data/align/test_aligned_sequences.fasta"), "missing")]
    Traceback (most recent call last):
      ...
    augur.align.AlignmentError: ERROR: reference missing not found in alignment
    '''
    seqs = {s.name:s for s in aln}
    if reference in seqs:
        ref_array = np.array(seqs[reference])
        if "-" not in ref_array:
            print("No gaps in alignment to trim (with respect to the reference, %s)"%reference)
        ungapped = ref_array!='-'
        ref_aln_array = np.array(aln)[:,ungapped]
    else:
        raise AlignmentError("ERROR: reference %s not found in alignment"%reference)

    if False in ungapped and insertion_csv:
        analyse_insertions(aln, ungapped, insertion_csv)

    out_seqs = []
    for seq, seq_array in zip(aln, ref_aln_array):
        seq.seq = Seq.Seq(''.join(seq_array))
        out_seqs.append(seq)

    if "-" in ref_array:
        print("Trimmed gaps in", reference, "from the alignment")

    return out_seqs

def analyse_insertions(aln, ungapped, insertion_csv):
    ## Gather groups (runs) of insertions:
    insertion_coords = [] # python syntax - e.g. [0, 3, 5] means indexes 0,1 & 2 are insertions (w.r.t. ref), to the right of 0-based ref pos 5
    _open_idx = False
    _ref_idx = -1
    for i, in_ref in enumerate(ungapped):
        if not in_ref and _open_idx is False:
            _open_idx = i # insertion run start
        elif in_ref and _open_idx is not False:
            insertion_coords.append([_open_idx, i, _ref_idx])# insertion run has finished
            _open_idx = False
        if in_ref:
            _ref_idx += 1
    if _open_idx is not False:
        insertion_coords.append([_open_idx, len(ungapped), _ref_idx])

    # For each run of insertions (w.r.t. reference) collect the insertions we have
    insertions = [defaultdict(list) for ins in insertion_coords]
    for idx, insertion_coord in enumerate(insertion_coords):
        for seq in aln:
            try:
                # biopython>=1.80 does not have seq.ungap()
                s = (seq[insertion_coord[0]:insertion_coord[1]].seq
                    .replace("-", "")
                    .replace("N", "")
                    .replace("?", "")
                )
            except AttributeError:
                # biopython<1.79 does not have seq.replace()
                s = (seq[insertion_coord[0]:insertion_coord[1]].seq
                    .ungap("-")
                    .ungap("N")
                    .ungap("?")
                )
            if len(s):
                insertions[idx][str(s)].append(seq.name)

    # output for auspice drag&drop -- GFF is 1-based & insertions are to the right of the base.
    header = ["strain"]+["insertion: {}bp @ ref pos {}".format(ic[1]-ic[0], ic[2]+1) for ic in insertion_coords]
    strain_data = defaultdict(lambda: ["" for _ in range(0, len(insertion_coords))])
    for idx, i_data in enumerate(insertions):
        for insertion_seq, strains in i_data.items():
            for strain in strains:
                strain_data[strain][idx] = insertion_seq
        if not len(i_data.keys()):
            # This happens when there _is_ an insertion, but it's an insertion of gaps
            # We know that these are associated with poor alignments...
            # GFF is 1-based & insertions are to the right of the base.
            insertion_coord = insertion_coords[idx]
            print("WARNING: {}bp insertion at ref position {} was due to 'N's or '?'s in provided sequences"\
                    .format(insertion_coord[1]-insertion_coord[0], insertion_coord[2]+1))
    with open_file(insertion_csv, 'w') as fh:
        print(",".join(header), file=fh)
        for strain in strain_data:
            print("{},{}".format(strain, ",".join(strain_data[strain])), file=fh)


def prettify_alignment(aln):
    '''
    Converts all bases to uppercase and removes auto reverse-complement prefix (_R_).
    This modifies the alignment in place.

    Parameters
    ----------
    aln : Bio.Align.MultipleSeqAlignment
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
    aln : Bio.Align.MultipleSeqAlignment
        Biopython Alignment
    '''
    for seq in aln:
        _seq = str(seq.seq)
        _seq = _seq.replace('-', 'N')
        seq.seq = Seq.Seq(_seq)


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
        elif isinstance(sample, (list, Align.MultipleSeqAlignment)):
            for s in sample:
                add(s.name)
        elif isinstance(sample, str):
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
    Return a set of seqs excluding those already in the alignment & print a warning
    message for each sequence which is exluded.
    """
    ret = []
    aln_names = {s.name for s in aln}
    for seq in seqs:
        if seq.name in aln_names:
            print("Excluding {} as it is already present in the alignment".format(seq.name))
        else:
            ret.append(seq)
    return ret
