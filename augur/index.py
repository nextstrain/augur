"""Count occurrence of bases/residues in a set of sequences.
"""

from itertools import combinations
import csv

from .io.file import open_file
from .io.sequences import read_sequences, is_vcf
from .errors import AugurError
from treetime.vcf_utils import read_vcf

DELIMITER = '\t'
ID_COLUMN = 'strain'

NUCLEOTIDE_CHARACTERS = {
    'A': 'a',
    'C': 'c',
    'G': 'g',
    'T': 't',
    'N': 'n',
    'other_IUPAC': 'ryswkmbdhv',
    '-': '-',
    '?': '?',
}

AMINO_ACID_CHARACTERS = {
    'A': 'a',  # Ala, Alanine
    'C': 'c',  # Cys, Cysteine
    'D': 'd',  # Asp, Aspartic acid
    'E': 'e',  # Glu, Glutamic acid
    'F': 'f',  # Phe, Phenylalanine
    'G': 'g',  # Gly, Glycine
    'H': 'h',  # His, Histidine
    'I': 'i',  # Ile, Isoleucine
    'K': 'k',  # Lys, Lysine
    'L': 'l',  # Leu, Leucine
    'M': 'm',  # Met, Methionine
    'N': 'n',  # Asn, Asparagine
    'P': 'p',  # Pro, Proline
    'Q': 'q',  # Gln, Glutamine
    'R': 'r',  # Arg, Arginine
    'S': 's',  # Ser, Serine
    'T': 't',  # Thr, Threonine
    'V': 'v',  # Val, Valine
    'W': 'w',  # Trp, Tryptophan
    'Y': 'y',  # Tyr, Tyrosine
    '*': '*',  # Stop codon
    'X': 'x',  # Unknown/ambiguous
    '-': '-',  # Gap
}



def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("index", help=__doc__)
    parser.add_argument('--sequences', '-s', required=True, help="sequences in FASTA or VCF formats. Augur will summarize the content of FASTA sequences and only report the names of strains found in a given VCF.")
    parser.add_argument('--output', '-o', help="tab-delimited file containing the number of bases per sequence in the given file. Output columns include strain, length, counts for each valid character, and a count of invalid characters. For nucleotide sequences: A, C, G, T, N, other_IUPAC, -, ?. For amino acid sequences: the 20 standard amino acids, *, X, -.", required=True)
    parser.add_argument('--seq-type', default='nuc', choices=['nuc', 'aa'], help="Sequence type: 'nuc' or 'aa'")
    parser.add_argument('--verbose', '-v', action="store_true", help="print index statistics to stdout")
    return parser


def index_vcf(vcf_path, index_path):
    """Create an index with a list of strain names from a given VCF. We do not
    calculate any statistics for VCFs.

    Parameters
    ----------
    vcf_path : str or `os.PathLike`
        path to a VCF file to index.
    index_path : str or `os.PathLike`
        path to a tab-delimited file containing the composition details for each
        sequence in the given input file.

    Returns
    -------
    int :
        number of strains indexed

    """
    strains = list(read_vcf(vcf_path)['sequences'].keys())

    num_of_seqs = 0

    with open_file(index_path, 'wt', newline='') as out_file:
        tsv_writer = csv.writer(out_file, delimiter = DELIMITER)

        #write header i output file
        header = [ID_COLUMN]
        tsv_writer.writerow(header)

        for record in strains:
            tsv_writer.writerow([record])
            num_of_seqs += 1

    return num_of_seqs


def index_sequence(sequence, values):
    """Count the number of characters for a given sequence record.

    Parameters
    ----------
    sequence : Bio.SeqRecord.SeqRecord
        sequence record to index.

    values : list of set of str
        values to count; sets must be non-overlapping and contain only
        single-character, lowercase strings

    Returns
    -------
    list :
        summary of the given sequence's strain name, length, character counts
        for the given values, and a final column with the number of characters
        that didn't match any of those in the given values.

    Examples
    --------
    >>> import Bio
    >>> nuc_values = [set(v) for v in NUCLEOTIDE_CHARACTERS.values()]
    >>> sequence_a = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq("ACTGN-?XWN"), id="seq_A")
    >>> index_sequence(sequence_a, nuc_values)
    ['seq_A', 10, 1, 1, 1, 1, 2, 1, 1, 1, 1]

    >>> sequence_b = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq("ACTGACTG"), id="seq_B")
    >>> index_sequence(sequence_b, nuc_values)
    ['seq_B', 8, 2, 2, 2, 2, 0, 0, 0, 0, 0]

    Characters in the given sequence that are not in the given list of values to
    count get counted in the final column.

    >>> sequence_c = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq("ACTG%@!!!NN"), id="seq_C")
    >>> index_sequence(sequence_c, nuc_values)
    ['seq_C', 11, 1, 1, 1, 1, 2, 0, 0, 0, 5]

    Amino acid sequences work the same way with amino acid character sets.

    >>> aa_values = [set(v) for v in AMINO_ACID_CHARACTERS.values()]
    >>> sequence_e = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq("MFLEKVG"), id="seq_E")
    >>> index_sequence(sequence_e, aa_values)
    ['seq_E', 7, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]

    Invalid amino acid characters are counted in the final column.

    >>> sequence_f = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq("MFLEK-X*?!"), id="seq_F")
    >>> index_sequence(sequence_f, aa_values)
    ['seq_F', 10, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2]

    The list of value sets must not overlap.

    >>> sequence_d = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq("A!C!TGXN"), id="seq_D")
    >>> index_sequence(sequence_d, [set('actg'), set('xn'), set('n')]) # doctest: +ELLIPSIS
    Traceback (most recent call last):
      ...
    ValueError: character sets ... and {'n'} overlap: {'n'}

    Value sets must contain only single-character, lowercase strings.

    >>> index_sequence(sequence_d, [{'a'}, {'c'}, {'T'}, {'g'}])
    Traceback (most recent call last):
      ...
    ValueError: character set {'T'} contains a non-lowercase character: 'T'

    >>> index_sequence(sequence_d, [{'actg'}])
    Traceback (most recent call last):
      ...
    ValueError: character set {'actg'} contains a multi-character (or maybe zero-length) string: 'actg'

    >>> index_sequence(sequence_d, [{'a', 'c'}, {0, 1}])
    Traceback (most recent call last):
      ...
    ValueError: character set {0, 1} contains a non-string element: 0
    """
    # Sets must be non-overlapping and contain only single-character, lowercase
    # strings, otherwise our assumptions are broken.  This is not much hardship
    # on the caller.
    for v in values:
        for c in v:
            if not isinstance(c, str):
                raise ValueError(f"character set {v!r} contains a non-string element: {c!r}")
            if len(c) != 1:
                raise ValueError(f"character set {v!r} contains a multi-character (or maybe zero-length) string: {c!r}")
            if c != c.lower():
                raise ValueError(f"character set {v!r} contains a non-lowercase character: {c!r}")

    for a, b in combinations(values, 2):
        if not a.isdisjoint(b):
            raise ValueError(f"character sets {a!r} and {b!r} overlap: {a & b!r}")

    counts = []
    seq = sequence.seq.lower()
    l = len(seq)

    for v in values:
        counts.append(sum(map(lambda x: seq.count(x), v)))

    invalid_characters = l-sum(counts)
    row = [sequence.id, l]+counts+[invalid_characters]

    return row


def get_characters_by_sequence_type(sequence_type):
    """Return character sets and labels for the given sequence type.

    Parameters
    ----------
    sequence_type : str
        ``'nuc'`` for nucleotide or ``'aa'`` for amino acid.

    Returns
    -------
    tuple :
        (values, labels) where values is a list of character sets and labels is
        a list of corresponding column names.
    """
    if sequence_type == 'aa':
        chars = AMINO_ACID_CHARACTERS
    else:
        chars = NUCLEOTIDE_CHARACTERS
    values = [set(v) for v in chars.values()]
    labels = list(chars.keys())
    return values, labels


def index_sequences(sequences_path, sequence_index_path, sequence_type):
    """Count the number of each valid character and invalid characters in a set
    of sequences and write the composition as a data frame to the given sequence
    index path.

    The sequence type (nucleotide or amino acid) is auto-detected from the first
    sequence. Nucleotide sequences count A, C, G, T, N, other IUPAC characters,
    gaps, and ambiguous characters. Amino acid sequences count the 20 standard
    amino acids, stop codons, X, and gaps.

    Parameters
    ----------
    sequences_path : str or `os.PathLike`
        path to a sequence file to index.

    sequence_index_path : str or `os.PathLike`
        path to a tab-delimited file containing the composition details for each
        sequence in the given input file.

    Returns
    -------
    int :
        number of sequences indexed

    int :
        total length of sequences indexed

    """
    seqs = read_sequences(sequences_path)

    tot_length = 0
    num_of_seqs = 0
    char_values, char_labels = get_characters_by_sequence_type(sequence_type)
    invalid = "invalid_residues" if sequence_type=="aa" else "invalid_nucleotides"
    header = [ID_COLUMN, 'length'] + char_labels + [invalid]

    with open_file(sequence_index_path, 'wt', newline='') as out_file:
        tsv_writer = csv.writer(out_file, delimiter = '\t', lineterminator='\n')
        tsv_writer.writerow(header)

        for record in seqs:
            row = index_sequence(record, char_values)
            tsv_writer.writerow(row)

            tot_length += row[1]
            num_of_seqs += 1

    return num_of_seqs, tot_length


def run(args):
    '''
    Index a set of sequences, counting character composition per sequence.
    Supports both nucleotide and amino acid sequences (auto-detected).
    '''
    if is_vcf(args.sequences):
        num_of_seqs = index_vcf(args.sequences, args.output)
        tot_length = None
    else:
        num_of_seqs, tot_length = index_sequences(args.sequences, args.output, args.seq_type)

    # empty input is fatal
    if not num_of_seqs:
        raise AugurError("Empty sequences provided")

    if args.verbose:
        if tot_length:
            print("Analysed %i sequences with an average length of %i characters." % (num_of_seqs, int(tot_length / num_of_seqs)))
        else:
            print("Analysed %i sequences" % num_of_seqs)
