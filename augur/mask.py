"""
Mask specified sites from a VCF or FASTA file.
"""
import os
import sys
from shutil import copyfile

import numpy as np
import pandas as pd
from Bio import SeqIO

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
            sites_to_mask.extend(list(range(int(row[1]), int(row[2])+1)))
        except ValueError as err:
            # Skip unparseable lines, including header lines.
            print("Could not read line %d of BED file %s: %s. Continuing." % (idx, mask_file, err))
    sites_to_mask = np.unique(sites_to_mask).tolist()
    print("Found %d sites to mask" % len(sites_to_mask))
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

    exclude = [chrom_name + "\t" + str(pos) for pos in mask_sites]
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

def mask_fasta(mask_sites, in_file, out_file):
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
    """
    # Load alignment as FASTA generator to prevent loading the whole alignment
    # into memory.
    alignment = SeqIO.parse(in_file, "fasta")

    # Write the masked alignment to disk one record at a time.
    print("Removing masked sites from FASTA file.")
    with open_file(out_file, "w") as oh:
        for record in alignment:
            # Convert to a mutable sequence to enable masking with Ns.
            sequence = record.seq.tomutable()
            sequence_length = len(sequence)
            # Replace all excluded sites with Ns.
            for site in mask_sites:
                if site < sequence_length:
                    sequence[site] = "N"
            record.seq = sequence
            SeqIO.write(record, oh, "fasta")

def register_arguments(parser):
    parser.add_argument('--sequences', '-s', required=True, help="sequences in VCF or FASTA format")
    parser.add_argument('--mask', required=True, help="locations to be masked in BED file format")
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
        return 1
    if not os.path.isfile(args.mask):
        print("ERROR: File {} does not exist!".format(args.mask))
        return 1
    if os.path.getsize(args.sequences) == 0:
        print("ERROR: {} is empty. Please check how this file was produced. "
              "Did an error occur in an earlier step?".format(args.sequences))
        return 1
    if os.path.getsize(args.mask) == 0:
        print("ERROR: {} is an empty file.".format(args.mask))
        return 1

    mask_sites = read_bed_file(args.mask)

    # For both FASTA and VCF masking, we need a proper separate output file
    if args.output is not None:
        out_file = args.output
    else:
        out_file = os.path.join(os.path.dirname(args.sequences),
                                "masked_" + os.path.basename(args.sequences))

    if is_vcf(args.sequences):
        mask_vcf(mask_sites, args.sequences, out_file, args.cleanup)
    else:
        mask_fasta(mask_sites, args.sequences, out_file)

    if args.output is None:
        copyfile(out_file, args.sequences)
        if args.cleanup:
            os.remove(out_file)
