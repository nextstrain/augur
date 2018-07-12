import os
from shutil import copyfile
import numpy as np
from Bio import AlignIO, SeqIO, Seq

def make_gaps_ambiguous(aln):
    '''
    replace all gaps by 'N' in all sequences in the alignment. TreeTime will treat them
    as fully ambiguous and replace then with the most likely state
    '''
    for seq in aln:
        seq_array = np.array(seq)
        gaps = seq_array == '-'
        seq_array[gaps] = 'N'
        seq.seq = Seq.Seq("".join(seq_array))


def run(args):
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
        print("Cannot read sequences -- make sure the file %s exists and contains sequences in"
              "fasta format"%seq_fname)
        return -1

    if ref_name and (ref_name not in seqs):
        print("Specified reference name %s is not in the sequence sample. Will not trim."%ref_name)
        ref_name = None

    # potentially add the reference seuqnce to the sequences
    if ref_fname and (not ref_name):
        if os.path.isfile(ref_fname):
            try:
                file_type = 'genbank' if ref_fname.split('.')[-1] in ['gb', 'genbank'] else 'fasta'
                ref_seq = SeqIO.read(ref_fname, file_type)
            except:
                print("WARNING: Cannot read reference sequence."
                      "\n\tMake sure the file %s contains one sequence in genbank or"
                      "fasta format"%ref_fname)
            else:
                ref_name = ref_seq.id
                seq_fname += '.ref.fasta'
                SeqIO.write(list(seqs.values())+[ref_seq], seq_fname, 'fasta')
        else:
            print("WARNING: Cannot read reference sequence."
                  "\n\tmake sure the file %s does not exist"%ref_fname)

    # before aligning, make a copy of the data that the aligner receives as input
    # (very useful for debugging purposes)
    copyfile(seq_fname, output+".pre_aligner.fasta")

    # align
    if args.method == 'mafft':
        cmd = """
              mafft --reorder --anysymbol --thread %d %s 1> %s 2> %s.log
              """%(args.nthreads, seq_fname, output, output)
        os.system(cmd)
        print("\nusing mafft to align via:\n\t" + cmd +
              " \n\n\tKatoh et al, Nucleic Acid Research, vol 30, issue 14"
              "\n\thttps://doi.org/10.1093%2Fnar%2Fgkf436\n")
    else:
        print('ERROR: alignment method not implemented')
        return -1

    # after aligning, make a copy of the data that the aligner produced (useful for debugging)
    copyfile(output, output+".post_aligner.fasta")

    # convert the aligner output to upper case (replacing the file in place)
    aln = AlignIO.read(output, 'fasta')
    for seq in aln:
        seq.seq = seq.seq.upper()
    AlignIO.write(aln, output, 'fasta')

    # if there's a valid reference in the alignment, strip out all the columns not present
    # in the reference, this will overwrite the alignment file
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
    '''
    aln = AlignIO.read(alignment_fname, 'fasta')
    seqs = {s.name:s for s in aln}
    if reference in seqs:
        ref_array = np.array(seqs[reference])
        ungapped = ref_array != '-'
        ref_aln_array = np.array(aln)[:, ungapped]
    else:
        print("WARNING: reference", reference, "not found in alignment")
        return

    out_seqs = []
    for seq, seq_array in zip(aln, ref_aln_array):
        seq.seq = Seq.Seq(''.join(seq_array))
        if keep_reference or seq.name != reference:
            out_seqs.append(seq)

    print("Trimmed gaps in", reference, "from the alignment")
    return out_seqs
