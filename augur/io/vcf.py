import os
import shlex

from .file import open_file
from .shell_command_runner import run_shell_command


shquote = shlex.quote


def is_vcf(filename):
    """Convenience method to check if a file is a vcf file.

    Examples
    --------
    >>> is_vcf(None)
    False
    >>> is_vcf("./foo")
    False
    >>> is_vcf("./foo.vcf")
    True
    >>> is_vcf("./foo.vcf.GZ")
    True
    """
    return bool(filename) and any(filename.lower().endswith(x) for x in ('.vcf', '.vcf.gz'))


def write_vcf(input_filename, output_filename, dropped_samps):
    if _filename_gz(input_filename):
        input_arg = "--gzvcf"
    else:
        input_arg = "--vcf"

    if _filename_gz(output_filename):
        output_pipe = "| gzip -c"
    else:
        output_pipe = ""

    drop_args = ["--remove-indv " + shquote(s) for s in dropped_samps]

    call = ["vcftools"] + drop_args + [input_arg, shquote(input_filename), "--recode --stdout", output_pipe, ">", shquote(output_filename)]

    print("Filtering samples using VCFTools with the call:")
    print(" ".join(call))
    run_shell_command(" ".join(call), raise_errors = True)
    # remove vcftools log file
    try:
        os.remove('out.log')
    except OSError:
        pass


def write_VCF_translation(prot_dict, vcf_file_name, ref_file_name):
    """
    Writes out a VCF-style file (which seems to be minimally handleable
    by vcftools and pyvcf) of the AA differences between sequences and the reference.
    This is a similar format created/used by read_in_vcf except that there is one
    of these dicts (with sequences, reference, positions) for EACH gene.

    Also writes out a fasta of the reference alignment.

    EBH 12 Dec 2017
    """
    import numpy as np

    #for the header
    seqNames = list(prot_dict[list(prot_dict.keys())[0]]['sequences'].keys())

    #prepare the header of the VCF & write out
    header=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+seqNames
    with open_file(vcf_file_name, 'w') as the_file:
        the_file.write( "##fileformat=VCFv4.2\n"+
                        "##source=NextStrain_Protein_Translation\n"+
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        the_file.write("\t".join(header)+"\n")

    refWrite = []
    vcfWrite = []

    #go through for every gene/protein
    for fname, prot in prot_dict.items():
        sequences = prot['sequences']
        ref = prot['reference']
        positions = prot['positions']

        #write out the reference fasta
        refWrite.append(">"+fname)
        refWrite.append(ref)

        #go through every variable position
        #There are no deletions here, so it's simpler than for VCF nuc sequenes!
        for pi in positions:
            pos = pi+1 #change numbering to match VCF not python
            refb = ref[pi] #reference base at this position

            #try/except is (much) faster than list comprehension!
            pattern = []
            for k,v in sequences.items():
                try:
                    pattern.append(sequences[k][pi])
                except KeyError:
                    pattern.append('.')
            pattern = np.array(pattern)

            #get the list of ALTs - minus any '.'!
            uniques = np.unique(pattern)
            uniques = uniques[np.where(uniques!='.')]

            #Convert bases to the number that matches the ALT
            j=1
            for u in uniques:
                pattern[np.where(pattern==u)[0]] = str(j)
                j+=1
            #Now convert these calls to #/# (VCF format)
            calls = [ j+"/"+j if j!='.' else '.' for j in pattern ]
            if len(uniques)==0:
                print("UNEXPECTED ERROR WHILE CONVERTING TO VCF AT POSITION {}".format(str(pi)))
                break

            #put it all together and write it out
            output = [fname, str(pos), ".", refb, ",".join(uniques), ".", "PASS", ".", "GT"] + calls

            vcfWrite.append("\t".join(output))

    #write it all out
    with open_file(ref_file_name, 'w') as the_file:
        the_file.write("\n".join(refWrite))

    with open_file(vcf_file_name, 'a') as the_file:
        the_file.write("\n".join(vcfWrite))

    if vcf_file_name.lower().endswith('.gz'):
        import os
        #must temporarily remove .gz ending, or gzip won't zip it!
        os.rename(vcf_file_name, vcf_file_name[:-3])
        call = ["gzip", vcf_file_name[:-3]]
        run_shell_command(" ".join(call), raise_errors = True)


def _filename_gz(filename):
    return filename.lower().endswith(".gz")
