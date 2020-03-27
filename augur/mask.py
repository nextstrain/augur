"""
Mask specified sites from a VCF file.
"""

import gzip
import os
from shutil import copyfile

import numpy as np
import pandas as pd
from Bio import SeqIO

from .utils import run_shell_command, shquote, open_file, is_vcf

def get_chrom_name(vcf_file):
    with open_file(vcf_file, mode='r') as f:
        for line in f:
            if line[0] != "#":
                header = line.strip().partition('\t')
                return header[0]

    print("ERROR: Something went wrong reading your VCF file: a CHROM column could not be found. "
          "Please check the file is valid VCF format.")
    return None


def read_bed_file(mask_file):
    #Read in BED file - 2nd column always chromStart, 3rd always chromEnd
    #I timed this against sets/update/sorted; this is faster
    sitesToMask = []
    bed = pd.read_csv(mask_file, sep='\t')
    for _, row in bed.iterrows():
        sitesToMask.extend(list(range(row[1], row[2]+1)))
    sitesToMask = np.unique(sitesToMask)
    return sitesToMask

def register_arguments(parser):
    parser.add_argument('--sequences', '-s', required=True, help="sequences in VCF format")
    parser.add_argument('--mask', required=True, help="locations to be masked in BED file format")
    parser.add_argument('--output', '-o', help="output file")

def mask_vcf(mask_sites, in_file, out_file, cleanup=True):
    cleanup_files = ['out.log']

    # Create the temporary masking file to be passed to VCFTools
    # Need CHROM name from VCF file:
    chrom_name = get_chrom_name(in_file)
    if chrom_name is None:
        return 1
    exclude = [chrom_name + "\t" + str(pos) for pos in mask_sites]
    temp_mask_file = in_file +"_maskTemp"
    with open_file(temp_mask_file, 'w') as fh:
        fh.write("\n".join(exclude))
    cleanup_files.append(temp_mask_file)

    #vcftools doesn't like input/output being the same file.
    #If no output specified, they will be, so use copy of input we'll delete later
    if out_file is None:
        out_file = in_file
        in_file = in_file + "_temp"
        copyfile(out_file, in_file)
        cleanup_files.append(in_file)

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


def run(args):
    '''
    mask specified sites from the VCF.
    this occurs by removing them entirely from the VCF, essentially making
    them identical to the reference at the locations

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

    if is_vcf(args.sequences):
        mask_vcf(mask_sites, args.sequences, args.output)
