"""
Mask specified sites from a VCF or FASTA file.
"""
import os
import sys
from shlex import quote as shquote
from shutil import copyfile

from Bio.Seq import MutableSeq

from .argparse_ import ExtendOverwriteDefault
from .io.file import open_file
from .io.sequences import read_sequences, write_sequences, is_vcf
from .io.shell_command_runner import run_shell_command
from .utils import load_mask_sites, VALID_NUCLEOTIDES

def get_chrom_name(vcf_file):
    """Read the CHROM field from the first non-header line of a vcf file.

    Returns:
    str or None: Either the CHROM field or None if no non-comment line could be found.
    """
    with open_file(vcf_file, mode='r') as f:
        for line in f:
            if line[0] != "#":
                header = line.strip().partition('\t')
                return header[0]

def mask_vcf(mask_sites, in_file, out_file, cleanup=True):
    """Mask the provided site list from a VCF file and write to a new file.

    This function relies on 'vcftools --exclude-positions' to mask the requested sites.

    Parameters
    ----------
    mask_sites: list of int
        A list of site indexes to exclude from the vcf.
    in_file: str
        The path to the vcf file you wish to mask.
    out_file: str
        The path to write the resulting vcf to
    cleanup: bool
        Clean up the intermediate files, including the VCFTools log and mask sites file
    """
    cleanup_files = ['out.log']

    # Create the temporary masking file to be passed to VCFTools
    # Need CHROM name from VCF file:
    chrom_name = get_chrom_name(in_file)
    if chrom_name is None:
        print("ERROR: Something went wrong reading your VCF file: a CHROM column could not be found. "
              "Please check the file is valid VCF format.")
        sys.exit(1)

    # mask_sites is zero-indexed, VCFTools expects 1-indexed.
    exclude = [chrom_name + "\t" + str(pos + 1) for pos in mask_sites]
    temp_mask_file = in_file + "_maskTemp"
    with open_file(temp_mask_file, 'w') as fh:
        fh.write("\n".join(exclude))
    cleanup_files.append(temp_mask_file)

    #Read in/write out according to file ending
    in_call = "--gzvcf" if in_file.lower().endswith(".gz") else "--vcf"
    out_call = "| gzip -c" if out_file.lower().endswith(".gz") else ""

    call = ["vcftools", "--exclude-positions", shquote(temp_mask_file), in_call, shquote(in_file), "--recode --stdout", out_call, ">", shquote(out_file)]
    print("Removing masked sites from VCF file using vcftools... this may take some time. Call:")
    print(" ".join(call))
    run_shell_command(" ".join(call), raise_errors = True)
    # remove vcftools log file
    if cleanup:
        for file in cleanup_files:
            try:
                os.remove(file)
            except OSError:
                pass


def mask_sequence(sequence, mask_sites, mask_from_beginning, mask_from_end, mask_invalid):
    """Mask characters at the given sites in a single sequence record, modifying the
    record in place.

    Parameters
    ----------
    sequence : Bio.SeqRecord.SeqRecord
        A sequence to be masked
    mask_sites: list of int
        A list of site indexes to exclude from the FASTA.
    mask_from_beginning: int
        Number of sites to mask from the beginning of each sequence (default 0)
    mask_from_end: int
        Number of sites to mask from the end of each sequence (default 0)
    mask_invalid: bool
        Mask invalid nucleotides (default False)

    Returns
    -------
    Bio.SeqRecord.SeqRecord
        Masked sequence in its original record object

    """
    # Convert to a mutable sequence to enable masking with Ns.
    sequence_length = len(sequence.seq)
    beginning, end = mask_from_beginning, mask_from_end

    if beginning + end > sequence_length:
        beginning, end = sequence_length, 0

    seq = str(sequence.seq)[beginning:-end or None]

    if mask_invalid:
        seq = "".join(nuc if nuc in VALID_NUCLEOTIDES else "N" for nuc in seq)

    masked_sequence = MutableSeq("N" * beginning + seq + "N" * end)

    # Replace all excluded sites with Ns.
    for site in mask_sites:
        if site < sequence_length:
            masked_sequence[site] = "N"

    sequence.seq = masked_sequence

    return sequence


def mask_fasta(mask_sites, in_file, out_file, mask_from_beginning=0, mask_from_end=0, mask_invalid=False):
    """Mask the provided site list from a FASTA file and write to a new file.

    Masked sites are overwritten as "N"s.

    Parameters
    ----------
    mask_sites: list of int
        A list of site indexes to exclude from the FASTA.
    in_file: str
        The path to the FASTA file you wish to mask.
    out_file: str
        The path to write the resulting FASTA to
    mask_from_beginning: int
       Number of sites to mask from the beginning of each sequence (default 0)
    mask_from_end: int
       Number of sites to mask from the end of each sequence (default 0)
    mask_invalid: bool
        Mask invalid nucleotides (default False)
    """
    # Load alignment as FASTA generator to prevent loading the whole alignment
    # into memory.
    alignment = read_sequences(in_file)

    # Write the masked alignment to disk one record at a time.
    print("Removing masked sites from FASTA file.")

    masked_sequences = (
        mask_sequence(
            sequence,
            mask_sites,
            mask_from_beginning,
            mask_from_end,
            mask_invalid,
        )
        for sequence in alignment
    )
    sequences_written = write_sequences(
        masked_sequences,
        out_file,
        "fasta"
    )

def register_arguments(parser):
    """
    Add arguments to parser.
    Kept as a separate function than `register_parser` to continue to support
    unit tests that use this function to create argparser.
    """
    parser.add_argument('--sequences', '-s', required=True, help="sequences in VCF or FASTA format")
    parser.add_argument('--mask', dest="mask_file", required=False, help="locations to be masked in either BED file format, DRM format, or one 1-indexed site per line.")
    parser.add_argument('--mask-from-beginning', type=int, default=0, help="FASTA Only: Number of sites to mask from beginning")
    parser.add_argument('--mask-from-end', type=int, default=0, help="FASTA Only: Number of sites to mask from end")
    parser.add_argument('--mask-invalid', action='store_true', help="FASTA Only: Mask invalid nucleotides")
    parser.add_argument("--mask-sites", nargs='+', action=ExtendOverwriteDefault, type = int,  help="1-indexed list of sites to mask")
    parser.add_argument('--output', '-o', help="output file")
    parser.add_argument('--no-cleanup', dest="cleanup", action="store_false",
                        help="Leave intermediate files around. May be useful for debugging")

def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("mask", help=__doc__)
    register_arguments(parser)
    return parser


def run(args):
    '''
    Mask specified sites from the VCF or FASTA.

    For VCF files, his occurs by removing them entirely from the VCF, essentially making
    them identical to the reference at the locations.

    For FASTA files, masked sites are replaced with "N".

    If users don't specify output, will overwrite the input file.
    '''
    # Check files exist and are not empty
    if not os.path.isfile(args.sequences):
        print("ERROR: File {} does not exist!".format(args.sequences))
        sys.exit(1)
    if os.path.getsize(args.sequences) == 0:
        print("ERROR: {} is empty. Please check how this file was produced. "
              "Did an error occur in an earlier step?".format(args.sequences))
        sys.exit(1)
    if args.mask_file:
        if not os.path.isfile(args.mask_file):
            print("ERROR: File {} does not exist!".format(args.mask_file))
            sys.exit(1)
        if os.path.getsize(args.mask_file) == 0:
            print("ERROR: {} is an empty file.".format(args.mask_file))
            sys.exit(1)
    if not any((args.mask_file, args.mask_from_beginning, args.mask_from_end, args.mask_sites, args.mask_invalid)):
        print("No masking sites provided. Must include one of --mask, --mask-from-beginning, --mask-from-end, --mask-invalid, or --mask-sites")
        sys.exit(1)

    mask_sites = set()
    if args.mask_sites:
        # Mask sites passed in as 1-indexed
        mask_sites.update(site - 1 for site in args.mask_sites)
    if args.mask_file:
        mask_sites.update(load_mask_sites(args.mask_file))
    mask_sites = sorted(mask_sites)

    # For both FASTA and VCF masking, we need a proper separate output file
    if args.output is not None:
        out_file = args.output
    else:
        out_file = os.path.join(os.path.dirname(args.sequences),
                                "masked_" + os.path.basename(args.sequences))

    if is_vcf(args.sequences):
        if args.mask_from_beginning or args.mask_from_end or args.mask_invalid:
            print("Cannot use --mask-from-beginning, --mask-from-end, or --mask-invalid with VCF files!")
            sys.exit(1)
        mask_vcf(mask_sites, args.sequences, out_file, args.cleanup)
    else:
        mask_fasta(mask_sites, args.sequences, out_file,
                   mask_from_beginning=args.mask_from_beginning,
                   mask_from_end=args.mask_from_end,
                   mask_invalid=args.mask_invalid)

    if args.output is None:
        copyfile(out_file, args.sequences)
        if args.cleanup:
            os.remove(out_file)
