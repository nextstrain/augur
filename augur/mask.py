"""
Mask specified sites from a VCF file.
"""

from Bio import SeqIO
import pandas as pd
import os
import numpy as np
from .utils import run_shell_command

def get_mask_sites(vcf_file, mask_file):
    '''
    Creates a temporary file in correct format for vcftools to use
    (two-column, tab-seperated: "chromName" "position")
    '''

    #Need CHROM name from VCF file:
    import gzip
    opn = gzip.open if vcf_file.lower().endswith('.gz') else open
    with opn(vcf_file, mode='rt') as f: #'rt' necessary for gzip
        for line in f:
            if line[0] != "#":
                header = line.strip().split('\t')
                chromName = header[0]
                break   # once chrom is found, no need to go through rest

    #Read in BED file - 2nd column always chromStart, 3rd always chromEnd
    #I timed this against sets/update/sorted; this is faster
    sitesToMask = []
    bed = pd.read_csv(mask_file, sep='\t')
    for index, row in bed.iterrows():
        sitesToMask.extend(list(range(row[1], row[2]+1)))
    sitesToMask = np.unique(sitesToMask)

    exclude = []
    for pos in sitesToMask:
        exclude.append(chromName+"\t"+str(pos))

    tempMaskFile = mask_file+"_maskTemp"
    with open(tempMaskFile, 'w') as the_file:
        the_file.write("\n".join(exclude))

    return tempMaskFile


def register_arguments(parser):
    parser.add_argument('--sequences', '-s', required=True, help="sequences in VCF format")
    parser.add_argument('--mask', required=True, help="locations to be masked in BED file format")
    parser.add_argument('--output', '-o', help="output file")


def run(args):
    '''
    mask specified sites from the VCF.
    this occurs by removing them entirely from the VCF, essentially making
    them identical to the reference at the locations

    If users don't specify output, will overwrite the input file.
    '''

    tempMaskFile = get_mask_sites(args.sequences, args.mask)

    #Read in/write out according to file ending
    inCall = "--gzvcf" if args.sequences.lower().endswith('.gz') else "--vcf"
    if args.output:
        outCall = "| gzip -c" if args.output.lower().endswith('.gz') else ""
    else:
        outCall = "| gzip -c" if args.sequences.lower().endswith('.gz') else ""

    #vcftools doesn't like input/output being the same file.
    #If no output specified, they will be, so use copy of input we'll delete later
    in_file = args.sequences
    out_file = args.output
    if not(args.output):
        from shutil import copyfile
        out_file = in_file
        in_file = args.sequences+"_temp"
        copyfile(args.sequences, in_file)

    call = ["vcftools", "--exclude-positions", tempMaskFile, inCall, in_file, "--recode --stdout", outCall, ">", out_file]
    print("Removing masked sites from VCF file using vcftools... this may take some time. Call:")
    print(" ".join(call))
    run_shell_command(" ".join(call), raise_errors = True)
    os.remove(tempMaskFile) #remove masking file
    # remove vcftools log file
    try:
        os.remove('out.log')
    except OSError:
        pass

    #remove copy of input if there was no output specified
    if not(args.output):
        os.remove(in_file)
