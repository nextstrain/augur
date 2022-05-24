#!/usr/bin/env python3
"""Interfaces for reading and writing data also known as input/output (I/O)
"""
import os
import shlex
import Bio.SeqIO
import Bio.SeqRecord
import sys
from contextlib import contextmanager
import pandas as pd
from pathlib import Path
from xopen import xopen

from augur.io_support.shell_command_runner import ShellCommandRunner


shquote = shlex.quote


@contextmanager
def open_file(path_or_buffer, mode="r", **kwargs):
    """Opens a given file path and returns the handle.

    Transparently handles compressed inputs and outputs.

    Parameters
    ----------
    path_or_buffer : str or Path-like or IO buffer
        Name of the file to open or an existing IO buffer

    mode : str
        Mode to open file (read or write)

    Returns
    -------
    IO
        File handle object

    """
    try:
        with xopen(path_or_buffer, mode, **kwargs) as handle:
            yield handle
    except TypeError:
        yield path_or_buffer


def read_metadata(metadata_file, id_columns=("strain", "name"), chunk_size=None):
    """Read metadata from a given filename and into a pandas `DataFrame` or
    `TextFileReader` object.

    Parameters
    ----------
    metadata_file : str
        Path to a metadata file to load.
    id_columns : list[str]
        List of possible id column names to check for, ordered by priority.
    chunk_size : int
        Size of chunks to stream from disk with an iterator instead of loading the entire input file into memory.

    Returns
    -------
    pandas.DataFrame or pandas.TextFileReader

    Raises
    ------
    KeyError :
        When the metadata file does not have any valid index columns.


    For standard use, request a metadata file and get a pandas DataFrame.

    >>> read_metadata("tests/functional/filter/metadata.tsv").index.values[0]
    'COL/FLR_00024/2015'

    Requesting an index column that doesn't exist should produce an error.

    >>> read_metadata("tests/functional/filter/metadata.tsv", id_columns=("Virus name",))
    Traceback (most recent call last):
      ...
    Exception: None of the possible id columns (('Virus name',)) were found in the metadata's columns ('strain', 'virus', 'accession', 'date', 'region', 'country', 'division', 'city', 'db', 'segment', 'authors', 'url', 'title', 'journal', 'paper_url')

    We also allow iterating through metadata in fixed chunk sizes.

    >>> for chunk in read_metadata("tests/functional/filter/metadata.tsv", chunk_size=5):
    ...     print(chunk.shape)
    ...
    (5, 14)
    (5, 14)
    (2, 14)

    """
    kwargs = {
        "sep": None,
        "engine": "python",
        "skipinitialspace": True,
        "na_filter": False,
    }

    if chunk_size:
        kwargs["chunksize"] = chunk_size

    # Inspect the first chunk of the metadata, to find any valid index columns.
    metadata = pd.read_csv(
        metadata_file,
        iterator=True,
        **kwargs,
    )
    chunk = metadata.read(nrows=1)
    metadata.close()

    id_columns_present = [
        id_column
        for id_column in id_columns
        if id_column in chunk.columns
    ]

    # If we couldn't find a valid index column in the metadata, alert the user.
    if not id_columns_present:
        raise Exception(f"None of the possible id columns ({id_columns!r}) were found in the metadata's columns {tuple(chunk.columns)!r}")
    else:
        index_col = id_columns_present[0]

    # If we found a valid column to index the DataFrame, specify that column and
    # also tell pandas that the column should be treated like a string instead
    # of having its type inferred. This latter argument allows users to provide
    # numerical ids that don't get converted to numbers by pandas.
    kwargs["index_col"] = index_col
    kwargs["dtype"] = {index_col: "string"}

    return pd.read_csv(
        metadata_file,
        **kwargs
    )


def read_sequences(*paths, format="fasta"):
    """Read sequences from one or more paths.

    Automatically infer compression mode (e.g., gzip, etc.) and return a stream
    of sequence records in the requested format (e.g., "fasta", "genbank", etc.).

    Parameters
    ----------
    paths : list of str or Path-like objects
        One or more paths to sequence files of any type supported by BioPython.

    format : str
        Format of input sequences matching any of those supported by BioPython
        (e.g., "fasta", "genbank", etc.).

    Yields
    ------
    Bio.SeqRecord.SeqRecord
        Sequence record from the given path(s).

    """
    for path in paths:
        # Open the given path as a handle, inferring the file's compression.
        # This way we can pass a handle to BioPython's SeqIO interface
        # regardless of the compression mode.
        with open_file(path) as handle:
            sequences = Bio.SeqIO.parse(handle, format)

            for sequence in sequences:
                yield sequence


def write_sequences(sequences, path_or_buffer, format="fasta"):
    """Write sequences to a given path in the given format.

    Automatically infer compression mode (e.g., gzip, etc.) based on the path's
    filename extension.

    Parameters
    ----------
    sequences : iterable of Bio.SeqRecord.SeqRecord objects
        A list-like collection of sequences to write

    path_or_buffer : str or Path-like object or IO buffer
        A path to a file to write the given sequences in the given format.

    format : str
        Format of input sequences matching any of those supported by BioPython
        (e.g., "fasta", "genbank", etc.)

    Returns
    -------
    int :
        Number of sequences written out to the given path.

    """
    with open_file(path_or_buffer, "wt") as handle:
        # Bio.SeqIO supports writing to the same handle multiple times for specific
        # file formats. For the formats we use, this function call should work for
        # both a newly opened file handle or one that is provided by the caller.
        # For more details see:
        # https://github.com/biopython/biopython/blob/25f5152f4aeefe184a323db25694fbfe0593f0e2/Bio/SeqIO/__init__.py#L233-L251
        sequences_written = Bio.SeqIO.write(
            sequences,
            handle,
            format
        )

    return sequences_written


def print_err(*args):
    print(*args, file=sys.stderr)


def is_vcf(filename):
    """Convenience method to check if a file is a vcf file.

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


def read_vcf(filename):
    if filename.lower().endswith(".gz"):
        import gzip
        file = gzip.open(filename, mode="rt", encoding='utf-8')
    else:
        file = open(filename, encoding='utf-8')

    chrom_line = next(line for line in file if line.startswith("#C"))
    file.close()
    headers = chrom_line.strip().split("\t")
    sequences = headers[headers.index("FORMAT") + 1:]

    # because we need 'seqs to remove' for VCF
    return sequences, sequences.copy()


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
    with open(vcf_file_name, 'w', encoding='utf-8') as the_file:
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
    with open(ref_file_name, 'w', encoding='utf-8') as the_file:
        the_file.write("\n".join(refWrite))

    with open(vcf_file_name, 'a', encoding='utf-8') as the_file:
        the_file.write("\n".join(vcfWrite))

    if vcf_file_name.lower().endswith('.gz'):
        import os
        #must temporarily remove .gz ending, or gzip won't zip it!
        os.rename(vcf_file_name, vcf_file_name[:-3])
        call = ["gzip", vcf_file_name[:-3]]
        run_shell_command(" ".join(call), raise_errors = True)


def run_shell_command(cmd, raise_errors=False, extra_env=None):
    """
    Run the given command string via Bash with error checking.

    Returns True if the command exits normally.  Returns False if the command
    exits with failure and "raise_errors" is False (the default).  When
    "raise_errors" is True, exceptions are rethrown.

    If an *extra_env* mapping is passed, the provided keys and values are
    overlayed onto the default subprocess environment.
    """
    return ShellCommandRunner(cmd, raise_errors=raise_errors, extra_env=extra_env).run()


def _filename_gz(filename):
    return filename.lower().endswith(".gz")
