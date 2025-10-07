import Bio.SeqIO
import os

from augur.errors import AugurError
from augur.utils import augur
from importlib.metadata import version as installed_version
from packaging.version import Version
from shlex import quote as shquote
from shutil import which
from tempfile import NamedTemporaryFile
from textwrap import dedent
from typing import Iterator, Iterable, Union
from .file import open_file
from .print import _n, indented_list
from .shell_command_runner import run_shell_command


def get_biopython_format(augur_format: str) -> str:
    """Validate sequence file format and return the inferred Biopython format."""
    supported_formats = {"fasta", "genbank"}
    if augur_format not in supported_formats:
        raise AugurError(f"Sequence file format {augur_format!r} is invalid. Must be one of {', '.join(repr(f) for f in sorted(supported_formats))}.")

    # Supported formats are valid Biopython formats
    biopython_format = augur_format

    # Allow comments in FASTA format using fasta-pearson in later biopython versions
    if Version(installed_version("biopython")) >= Version("1.85") and biopython_format == "fasta":
        biopython_format = "fasta-pearson"

    return biopython_format


def read_sequences(
    *paths: Iterable[Union[str, os.PathLike]],
    format: str = "fasta",
) -> Iterator[Bio.SeqIO.SeqRecord]:
    """Read sequences from one or more paths.

    Automatically infer compression mode (e.g., gzip, etc.) and return a stream
    of sequence records given the file format.

    Parameters
    ----------
    paths
        One or more paths to sequence files.

    format
        Format of input sequences. Either "fasta" or "genbank".

    Returns
    -------
        Sequence records from the given path(s).
    """
    for path in paths:
        # Open the given path as a handle, inferring the file's compression.
        # This way we can pass a handle to BioPython's SeqIO interface
        # regardless of the compression mode.
        with open_file(path) as handle:
            sequences = Bio.SeqIO.parse(handle, get_biopython_format(format))

            for sequence in sequences:
                yield sequence


def read_single_sequence(
    path: Union[str, os.PathLike],
    format: str = "fasta",
) -> Bio.SeqIO.SeqRecord:
    """Read a single sequence from a path.

    Automatically infers compression mode.

    Parameters
    ----------
    path
        Path to a sequence file.

    format
        Format of input file. Either "fasta" or "genbank".

    Returns
    -------
        A single sequence record from the given path.
    """
    with open_file(path) as handle:
        return Bio.SeqIO.read(handle, get_biopython_format(format))


def write_sequences(sequences, path_or_buffer, format="fasta"):
    """Write sequences to a given path in the given format.

    Automatically infer compression mode (e.g., gzip, etc.) based on the path's
    filename extension.

    Parameters
    ----------
    sequences : iterable of Bio.SeqRecord.SeqRecord
        A list-like collection of sequences to write

    path_or_buffer : str or `os.PathLike` or `io.StringIO`
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


def write_records_to_fasta(records, fasta, seq_id_field='strain', seq_field='sequence'):
    """
    Write sequences from dict *records* to a *fasta* file.
    Yields the records with the *seq_field* dropped so that they can be consumed downstream.

    Parameters
    ----------
    records: iterable of dict
        Iterator that yields dict that contains sequences

    fasta: str
        Path to FASTA file

    seq_id_field: str, optional
        Field name for the sequence identifier

    seq_field: str, optional
        Field name for the genomic sequence

    Yields
    ------
    dict:
        A copy of the record with *seq_field* dropped

    Raises
    ------
    AugurError
        When the sequence id field or sequence field does not exist in a record
    """
    with open_file(fasta, "w") as output_fasta:
        for record in records:
            if seq_id_field not in record:
                raise AugurError(f"Provided sequence identifier field {seq_id_field!r} does not exist.")
            if seq_field not in record:
                raise AugurError(f"Provided sequence field {seq_field!r} does not exist.")

            output_fasta.writelines([
                f">{record[seq_id_field]}\n",
                f"{record[seq_field]}\n"
            ])

            yield {
                key: value
                for key, value in record.items()
                if key != seq_field
            }


def subset_fasta(input_filename: str, output_filename: str, ids_file: str, nthreads: int):
    command = f"""
        {seqkit()} --threads {nthreads} grep -f {ids_file} {shquote(input_filename)} |
        {augur()} write-file {shquote(output_filename)}
    """

    try:
        run_shell_command(command, raise_errors=True)
    except Exception:
        if os.path.isfile(output_filename):
            # Remove the partial output file.
            os.remove(output_filename)
            raise AugurError(f"Sequence output failed, see error(s) above.")
        else:
            raise AugurError(f"Sequence output failed, see error(s) above. The command may have already written data to stdout. You may want to clean up any partial outputs.")


def load_features(reference, feature_names=None):
    """
    Parse a GFF/GenBank reference file. See the docstrings for _read_gff and
    _read_genbank for details.

    Parameters
    ----------
    reference : str
        File path to GFF or GenBank (.gb) reference
    feature_names : None or set or list (optional)
        Restrict the genes we read to those in the set/list

    Returns
    -------
    features : dict
        keys: feature names, values: :py:class:`Bio.SeqFeature.SeqFeature` Note
        that feature names may not equivalent to GenBank feature keys

    Raises
    ------
    AugurError
        If the reference file doesn't exist, is malformed / empty, or has errors
    """
    #checks explicitly for GFF otherwise assumes Genbank
    if not os.path.isfile(reference):
        raise AugurError(f"reference sequence file {reference!r} not found")

    if '.gff' in reference.lower():
        features = _read_gff(reference, feature_names)
    else:
        features = _read_genbank(reference, feature_names)

    if errors := find_feature_errors(features):
        raise AugurError(dedent(f"""\
            Reference file {reference!r} has errors:
                {indented_list(errors, "                ")}"""))

    return features


def find_feature_errors(features):
    """
    Find and return errors for features parsed from a GFF/GenBank reference
    file.

    Parameters
    ----------
    features : dict
        keys: feature names, values: :py:class:`Bio.SeqFeature.SeqFeature`

    Returns
    -------
    list of str
        Error messages
    """
    errors = []

    for feature_name, feat in features.items():
        if feature_name == 'nuc':
            continue

        # Error if feature length is not a multiple of 3.
        length = len(feat.location)
        if length % 3:
            errors.append(f"{feature_name!r} has length {length} which is not a multiple of 3.")

    return errors


def _read_nuc_annotation_from_gff(record, reference):
    """
    Looks for the ##sequence-region pragma as well as 'region' & 'source' GFF
    types. Note that 'source' isn't really a GFF feature type, but is used
    widely in the Nextstrain ecosystem. If there are multiple we check that the
    coordinates agree.

    Parameters
    ----------
    record : :py:class:`Bio.SeqRecord.SeqRecord`
    reference: string
        File path to GFF reference

    Returns
    -------
    :py:class:`Bio.SeqFeature.SeqFeature`

    Raises
    ------
    AugurError
        If no information on the genome / seqid length is available or if the
        information is contradictory
    """
    nuc = {}
    # Attempt to parse the sequence-region pragma to learn the genome
    # length (in the absence of record/source we'll use this for 'nuc')
    sequence_regions = record.annotations.get('sequence-region', [])
    if len(sequence_regions)>1:
        raise AugurError(f"Reference {reference!r} contains multiple ##sequence-region pragma lines. Augur can only handle GFF files with a single one.")
    elif sequence_regions:
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        (name, start, stop) = sequence_regions[0]
        nuc['pragma'] = SeqFeature(
            FeatureLocation(start, stop, strand=1),
            type='##sequence-region pragma',
            id=name,
        )
    for feat in record.features:
        if feat.type == "region":
            nuc['region'] = feat
        elif feat.type == "source":
            nuc['source'] = feat

    # ensure they all agree on coordinates, if there are multiple
    if len(nuc.values())>1:
        coords = [(name, int(feat.location.start), int(feat.location.end)) for name,feat in nuc.items()]
        if not all(el[1]==coords[0][1] and el[2]==coords[0][2] for el in coords):
            raise AugurError(f"Reference {reference!r} contained contradictory coordinates for the seqid/genome. We parsed the following coordinates: " +
                             ', '.join([f"{el[0]}: [{el[1]+1}, {el[2]}]" for el in coords]) # +1 on the first coord to shift to one-based GFF representation
                             )

    if 'pragma' in nuc: ## the pragma is GFF's preferred way to define nuc coords
        return nuc['pragma']
    elif 'region' in nuc:
        return nuc['region']
    elif 'source' in nuc:
        return nuc['source']
    else:
        raise AugurError(f"Reference {reference!r} didn't define any information we can use to create the 'nuc' annotation. You can use a line with a 'record' or 'source' GFF type or a ##sequence-region pragma.")


def _read_gff(reference, feature_names):
    """
    Read a GFF file. We only read GFF IDs 'gene' or 'source' (the latter may not technically
    be a valid GFF field, but is used widely within the Nextstrain ecosystem).
    Only the first entry in the GFF file is parsed.
    We create a "feature name" via:
    - for 'source' IDs use 'nuc'
    - for 'gene' IDs use the 'gene', 'gene_name' or 'locus_tag'.
      If none are specified, the intention is to silently ignore but there are bugs here.

    Parameters
    ----------
    reference : string
        File path to GFF reference
    feature_names : None or set or list
        Restrict the genes we read to those in the set/list

    Returns
    -------
    features : dict
        keys: feature names, values: :py:class:`Bio.SeqFeature.SeqFeature`
        Note that feature names may not equivalent to GenBank feature keys

    Raises
    ------
    AugurError
        If the reference file contains no IDs or multiple different seqids
        If a gene is found with the name 'nuc'
    """
    from BCBio import GFF
    valid_types = ['gene', 'source', 'region']
    features = {}

    with open_file(reference) as in_handle:
        # Note that `GFF.parse` doesn't always yield GFF records in the order
        # one may expect, but since we raise AugurError if there are multiple
        # this doesn't matter.
        # TODO: Remove warning suppression after it's addressed upstream:
        # <https://github.com/chapmanb/bcbb/issues/143>
        import warnings
        from Bio import BiopythonDeprecationWarning
        warnings.simplefilter("ignore", BiopythonDeprecationWarning)
        gff_entries = list(GFF.parse(in_handle, limit_info={'gff_type': valid_types}))
        warnings.simplefilter("default", BiopythonDeprecationWarning)

        if len(gff_entries) == 0:
            raise AugurError(f"Reference {reference!r} contains no valid data rows. Valid GFF types (3rd column) are {', '.join(valid_types)}.")
        elif len(gff_entries) > 1:
            raise AugurError(f"Reference {reference!r} contains multiple seqids (first column). Augur can only handle GFF files with a single seqid.")
        else:
            rec = gff_entries[0]

        features['nuc'] = _read_nuc_annotation_from_gff(rec, reference)
        features_skipped = 0

        for feat in rec.features:
            if feat.type == "gene":
                # Check for gene names stored in qualifiers commonly used by
                # virus-specific gene maps first (e.g., 'gene',
                # 'gene_name'). Then, check for qualifiers used by non-viral
                # pathogens (e.g., 'locus_tag').
                if "gene" in feat.qualifiers:
                    fname = feat.qualifiers["gene"][0]
                elif "gene_name" in feat.qualifiers:
                    fname = feat.qualifiers["gene_name"][0]
                elif "locus_tag" in feat.qualifiers:
                    fname = feat.qualifiers["locus_tag"][0]
                else:
                    features_skipped+=1
                    fname = None

                if fname == 'nuc':
                    raise AugurError(f"Reference {reference!r} contains a gene with the name 'nuc'. This is not allowed.")

                if feature_names is not None and fname not in feature_names:
                    # Skip (don't store) this feature
                    continue

                if fname:
                    features[fname] = feat

        if feature_names is not None:
            for fe in feature_names:
                if fe not in features:
                    print("Couldn't find gene {} in GFF or GenBank file".format(fe))

        if features_skipped:
            print(f"WARNING: {features_skipped} GFF rows of type=gene skipped as they didn't have a gene, gene_name or locus_tag attribute.")

    return features


def _read_nuc_annotation_from_genbank(record, reference):
    """
    Extracts the mandatory 'source' feature. If the sequence is present we check
    the length agrees with the source. (The 'ORIGIN' may be left blank,
    according to <https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html>.)

    See <https://www.insdc.org/submitting-standards/feature-table/> for more.

    Parameters
    ----------
    record : :py:class:`Bio.SeqRecord.SeqRecord` reference: string

    Returns
    -------
    :py:class:`Bio.SeqFeature.SeqFeature`

    Raises
    ------
    AugurError
        If 'source' not defined or if coords contradict.
    """
    nuc = None
    for feat in record.features:
        if feat.type=='source':
            nuc = feat
    if not nuc:
        raise AugurError(f"Reference {reference!r} did not define the mandatory source feature.")
    if nuc.location.start!=0: # this is a '1' in the GenBank file
        raise AugurError(f"Reference {reference!r} source feature did not start at 1.")
    if record.seq and len(record.seq)!=nuc.location.end:
        raise AugurError(f"Reference {reference!r} source feature was length {nuc.location.end} but the included sequence was length {len(record.seq)}.")
    return nuc


def _read_genbank(reference, feature_names):
    """
    Read a GenBank file. We only read GenBank feature keys 'CDS' or 'source'.
    We create a "feature name" via:
    - for 'source' features use 'nuc'
    - for 'CDS' features use the locus_tag or the gene. If neither, then silently ignore.

    Parameters
    ----------
    reference : string
        File path to GenBank reference
    feature_names : None or set or list
        Restrict the CDSs we read to those in the set/list

    Returns
    -------
    features : dict
        keys: feature names, values: :py:class:`Bio.SeqFeature.SeqFeature`
        Note that feature names may not equivalent to GenBank feature keys

    Raises
    ------
    AugurError
        If 'nuc' annotation not parsed
        If a CDS feature is given the name 'nuc'
    """
    gb = read_single_sequence(reference, format='genbank')
    features = {
        'nuc': _read_nuc_annotation_from_genbank(gb, reference)
    }

    features_skipped = 0
    for feat in gb.features:
        if feat.type=='CDS':
            fname = None
            if "locus_tag" in feat.qualifiers:
                fname = feat.qualifiers["locus_tag"][0]
            elif "gene" in feat.qualifiers:
                fname = feat.qualifiers["gene"][0]
            else:
                features_skipped+=1

            if fname == 'nuc':
                raise AugurError(f"Reference {reference!r} contains a CDS with the name 'nuc'. This is not allowed.")

            if fname and (feature_names is None or fname in feature_names):
                features[fname] = feat

    if features_skipped:
        print(f"WARNING: {features_skipped} CDS features skipped as they didn't have a locus_tag or gene qualifier.")

    return features


# TODO: consider consolidating with augur.io.metadata.read_metadata_with_sequences
# <https://github.com/nextstrain/augur/pull/1821#discussion_r2138629823>
def read_sequence_ids(file: str, nthreads: int):
    """Get unique identifiers from a sequence file."""
    unique = set()
    duplicates = set()

    with NamedTemporaryFile("w+") as temp_file:
        if is_vcf(file):
            # Unfortunately, vcf-query doesn't support piped input from augur read-file.
            command = f"""
                vcf-query -l {shquote(file)} > {shquote(temp_file.name)}
            """
        else:
            command = f"""
                {seqkit()} --threads {nthreads} fx2tab --name --only-id {shquote(file)} > {shquote(temp_file.name)}
            """

        try:
            run_shell_command(command, raise_errors=True)
        except Exception as error:
            raise AugurError(f"Unable to read {file!r}. See error above.") from error

        temp_file.seek(0)

        for line in temp_file:
            identifier = line.strip()
            if identifier in unique:
                duplicates.add(identifier)
            else:
                unique.add(identifier)

    if duplicates:
        raise AugurError(dedent(f"""\
            Sequence ids must be unique.

            The following {_n("id was", "ids were", len(duplicates))} were duplicated in {file!r}:

              {indented_list(sorted(duplicates), '            ' + '  ')}
        """))
    
    return unique


def seqkit():
    """
    Internal helper for invoking SeqKit.

    Unlike ``augur.merge.sqlite3()``, this function is not a wrapper around
    subprocess.run. It is meant to be called without any parameters and only
    returns the location of the executable. This is due to differences in the
    way the two programs are invoked.
    """
    seqkit = os.environ.get("SEQKIT", which("seqkit"))
    if not seqkit:
        raise AugurError(dedent(f"""\
            Unable to find the program `seqkit`.  Is it installed?
            In order to handle FASTA files, SeqKit must be installed
            separately.  It is typically provided by a Nextstrain runtime.
            """))
    return shquote(seqkit)


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
