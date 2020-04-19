"""
Mask specified sites from a VCF or FASTA file.
"""
import os
import sys
from shutil import copyfile

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import MutableSeq

from .utils import run_shell_command, shquote, open_file, is_vcf

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

def read_bed_file(mask_file):
    """Read the full list of excluded sites from the BED file.

    Second column is chromStart, 3rd is chromEnd. Generate a range from these two columns.
    """
    sites_to_mask = []
    bed = pd.read_csv(mask_file, sep='\t', header=None, usecols=[1,2])
    for idx, row in bed.iterrows():
        try:
            sites_to_mask.extend(range(int(row[1]), int(row[2])))
        except ValueError as err:
            # Skip unparseable lines, including header lines.
            print("Could not read line %d of BED file %s: %s. Continuing." % (idx, mask_file, err))
    sites_to_mask = np.unique(sites_to_mask).tolist()
    print("Found %d sites to mask in '%s'" % (len(sites_to_mask), mask_file))
    return sites_to_mask

def mask_vcf(mask_sites, in_file, out_file, cleanup=True):
    """Mask the provided site list from a VCF file and write to a new file.

    This function relies on 'vcftools --exclude-positions' to mask the requested sites.

    Parameters:
    -----------
    mask_sites: list[int]
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

def mask_fasta(mask_sites, in_file, out_file, mask_from_beginning=0, mask_from_end=0):
    """Mask the provided site list from a FASTA file and write to a new file.

    Masked sites are overwritten as "N"s.

:
    -----------
    mask_sites: list[int]
        A list of site indexes to exclude from the FASTA.
    in_file: str
        The path to the FASTA file you wish to mask.
    out_file: str
        The path to write the resulting FASTA to
    mask_from_beginning: int
       Number of sites to mask from the beginning of each sequence (default 0)
    mask_from_end: int
       Number of sites to mask from the end of each sequence (default 0)
    """
    # Load alignment as FASTA generator to prevent loading the whole alignment
    # into memory.
    alignment = SeqIO.parse(in_file, "fasta")

    # Write the masked alignment to disk one record at a time.
    print("Removing masked sites from FASTA file.")
    with open_file(out_file, "w") as oh:
        for record in alignment:
            # Convert to a mutable sequence to enable masking with Ns.
            sequence_length = len(record.seq)
            beginning, end = mask_from_beginning, mask_from_end
            if beginning + end > sequence_length:
                beginning, end = sequence_length, 0
            sequence = MutableSeq(
                "N" * beginning +
                str(record.seq)[beginning:-end or None] +
                "N" * end
            )
            # Replace all excluded sites with Ns.
            for site in mask_sites:
                if site < sequence_length:
                    sequence[site] = "N"
            record.seq = sequence
            SeqIO.write(record, oh, "fasta")

def register_arguments(parser):
    parser.add_argument('--sequences', '-s', required=True, help="sequences in VCF or FASTA format")
    parser.add_argument('--mask', dest="mask_file", required=False, help="locations to be masked in BED file format")
    parser.add_argument('--mask-from-beginning', type=int, default=0, help="FASTA Only: Number of sites to mask from beginning")
    parser.add_argument('--mask-from-end', type=int, default=0, help="FASTA Only: Number of sites to mask from end")
    parser.add_argument("--mask-sites", nargs='+', type = int,  help="1-indexed list of sites to mask")
    parser.add_argument('--output', '-o', help="output file")
    parser.add_argument('--no-cleanup', dest="cleanup", action="store_false",
                        help="Leave intermediate files around. May be useful for debugging")

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
    if not any((args.mask_file, args.mask_from_beginning, args.mask_from_end, args.mask_sites)):
        print("No masking sites provided. Must include one of --mask, --mask-from-beginning, --mask-from-end, or --mask-sites")
        sys.exit(1)

    mask_sites = set()
    if args.mask_sites:
        # Mask sites passed in as 1-indexed
        mask_sites.update(site - 1 for site in args.mask_sites)
    if args.mask_file:
        mask_sites.update(read_bed_file(args.mask_file))
    mask_sites = sorted(mask_sites)

    # For both FASTA and VCF masking, we need a proper separate output file
    if args.output is not None:
        out_file = args.output
    else:
        out_file = os.path.join(os.path.dirname(args.sequences),
                                "masked_" + os.path.basename(args.sequences))

    if is_vcf(args.sequences):
        if args.mask_from_beginning or args.mask_from_end:
            print("Cannot use --mask-from-beginning or --mask-from-end with VCF files!")
            sys.exit(1)
        mask_vcf(mask_sites, args.sequences, out_file, args.cleanup)
    else:
        mask_fasta(mask_sites, args.sequences, out_file, 
                   mask_from_beginning=args.mask_from_beginning,
                   mask_from_end=args.mask_from_end)

    if args.output is None:
        copyfile(out_file, args.sequences)
        if args.cleanup:
            os.remove(out_file)
