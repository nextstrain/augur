import argparse
import Bio
import Bio.Phylo
from datetime import datetime
import gzip
import os, json, sys
import pandas as pd
import subprocess
import shlex
from contextlib import contextmanager
from treetime.utils import numeric_date
from collections import defaultdict
from pkg_resources import resource_stream
from io import TextIOWrapper
from .__version__ import __version__

from augur.io import open_file

from augur.util_support.color_parser import ColorParser
from augur.util_support.date_disambiguator import DateDisambiguator
from augur.util_support.metadata_file import MetadataFile
from augur.util_support.node_data_reader import NodeDataReader
from augur.util_support.shell_command_runner import ShellCommandRunner


class AugurException(Exception):
    pass


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

def myopen(fname, mode):
    if fname.endswith('.gz'):
        import gzip
        return gzip.open(fname, mode, encoding='utf-8')
    else:
        return open(fname, mode, encoding='utf-8')

def get_json_name(args, default=None):
    if args.output_node_data:
        return args.output_node_data
    else:
        if default:
            print("WARNING: no name for the output file was specified. Writing results to %s."%default, file=sys.stderr)
            return default
        else:
            raise ValueError("Please specify a name for the JSON file containing the results.")


def ambiguous_date_to_date_range(uncertain_date, fmt, min_max_year=None):
    return DateDisambiguator(uncertain_date, fmt=fmt, min_max_year=min_max_year).range()

def read_metadata(fname, query=None, as_data_frame=False):
    return MetadataFile(fname, query, as_data_frame).read()

def is_date_ambiguous(date, ambiguous_by="any"):
    """
    Returns whether a given date string in the format of YYYY-MM-DD is ambiguous by a given part of the date (e.g., day, month, year, or any parts).

    Parameters
    ----------
    date : str
        Date string in the format of YYYY-MM-DD
    ambiguous_by : str
        Field of the date string to test for ambiguity ("day", "month", "year", "any")
    """
    date_components = date.split('-', 2)

    if len(date_components) == 3:
        year, month, day = date_components
    elif len(date_components) == 2:
        year, month = date_components
        day = "XX"
    else:
        year = date_components[0]
        month = "XX"
        day = "XX"

    # Determine ambiguity hierarchically such that, for example, an ambiguous
    # month implicates an ambiguous day even when day information is available.
    return any((
        "X" in year,
        "X" in month and ambiguous_by in ("any", "month", "day"),
        "X" in day and ambiguous_by in ("any", "day")
    ))

def get_numerical_date_from_value(value, fmt=None, min_max_year=None, raise_error=True):
    if type(value)!=str:
        if raise_error:
            raise ValueError(value)
        else:
            numerical_date = None
    elif 'XX' in value:
        ambig_date = ambiguous_date_to_date_range(value, fmt, min_max_year)
        if ambig_date is None or None in ambig_date:
            numerical_date = [None, None] #don't send to numeric_date or will be set to today
        else:
            numerical_date = [numeric_date(d) for d in ambig_date]
    else:
        try:
            numerical_date = numeric_date(datetime.strptime(value, fmt))
        except:
            numerical_date = None

    return numerical_date

def get_numerical_dates(meta_dict, name_col = None, date_col='date', fmt=None, min_max_year=None):
    if fmt:
        numerical_dates = {}

        if isinstance(meta_dict, dict):
            for k,m in meta_dict.items():
                v = m[date_col]
                try:
                    numerical_dates[k] = get_numerical_date_from_value(
                        v,
                        fmt,
                        min_max_year
                    )
                except ValueError:
                    print(
                        "WARNING: %s has an invalid data string: %s"% (k, v),
                        file=sys.stderr
                    )
                    continue
        elif isinstance(meta_dict, pd.DataFrame):
            strains = meta_dict.index.values
            dates = meta_dict[date_col].apply(
                lambda date: get_numerical_date_from_value(
                    date,
                    fmt,
                    min_max_year,
                    raise_error=False
                )
            ).values
            numerical_dates = dict(zip(strains, dates))
    else:
        if isinstance(meta_dict, dict):
            numerical_dates = {k:float(v) for k,v in meta_dict.items()}
        elif isinstance(meta_dict, pd.DataFrame):
            strains = meta_dict.index.values
            dates = meta_dict[date_col].astype(float)
            numerical_dates = dict(zip(strains, dates))

    return numerical_dates


class InvalidTreeError(Exception):
    """Represents an error loading a phylogenetic tree from a filename.
    """
    pass


def read_tree(fname, min_terminals=3):
    """Safely load a tree from a given filename or raise an error if the file does
    not contain a valid tree.

    Parameters
    ----------
    fname : str
        name of a file containing a phylogenetic tree

    min_terminals : int
        minimum number of terminals required for the parsed tree as a sanity
        check on the tree

    Raises
    ------
    InvalidTreeError
        If the given file exists but does not seem to contain a valid tree format.

    Returns
    -------
    Bio.Phylo :
        BioPython tree instance

    """
    T = None
    supported_tree_formats = ["newick", "nexus"]
    for fmt in supported_tree_formats:
        try:
            T = Bio.Phylo.read(fname, fmt)

            # Check the sanity of the parsed tree to handle cases when non-tree
            # data are still successfully parsed by BioPython. Too few terminals
            # in a tree indicates that the input is not valid.
            if T.count_terminals() < min_terminals:
                T = None
            else:
                break
        except ValueError:
            # We cannot open the tree in the current format, so we will try
            # another.
            pass

    # If the tree cannot be loaded, raise an error to that effect.
    if T is None:
        raise InvalidTreeError(
            "Could not read the given tree %s using the following supported formats: %s" % (fname, ", ".join(supported_tree_formats))
        )

    return T


def read_node_data(fnames, tree=None):
    return NodeDataReader(fnames, tree).read()


def write_json(data, file_name, indent=(None if os.environ.get("AUGUR_MINIFY_JSON") else 2), include_version=True):
    """
    Write ``data`` as JSON to the given ``file_name``, creating parent directories
    if necessary. The augur version is included as a top-level key "augur_version".

    Parameters
    ----------
    data : dict
        data to write out to JSON
    file_name : str
        file name to write to
    indent : int or None, optional
        JSON indentation level. Default is `None` if the environment variable `AUGUR_MINIFY_JSON`
        is truthy, else 1
    include_version : bool, optional
        Include the augur version. Default: `True`.

    Raises
    ------
    OSError
    """
    #in case parent folder does not exist yet
    parent_directory = os.path.dirname(file_name)
    if parent_directory and not os.path.exists(parent_directory):
        try:
            os.makedirs(parent_directory)
        except OSError: #Guard against race condition
            if not os.path.isdir(parent_directory):
                raise

    if include_version:
        data["generated_by"] = {"program": "augur", "version": get_augur_version()}

    with open(file_name, 'w', encoding='utf-8') as handle:
        json.dump(data, handle, indent=indent, sort_keys=True)


def load_features(reference, feature_names=None):
    #read in appropriately whether GFF or Genbank
    #checks explicitly for GFF otherwise assumes Genbank
    if not os.path.isfile(reference):
        print("ERROR: reference sequence not found. looking for", reference)
        return None

    features = {}
    if '.gff' in reference.lower():
        #looks for 'gene' and 'gene' as best for TB
        try:
            from BCBio import GFF #Package name is confusing - tell user exactly what they need!
        except ImportError:
            print("ERROR: Package BCBio.GFF not found! Please install using \'pip install bcbio-gff\' before re-running.")
            return None
        limit_info = dict( gff_type = ['gene'] )

        with open(reference, encoding='utf-8') as in_handle:
            for rec in GFF.parse(in_handle, limit_info=limit_info):
                for feat in rec.features:
                    if feature_names is not None: #check both tags; user may have used either
                        if "gene" in feat.qualifiers and feat.qualifiers["gene"][0] in feature_names:
                            fname = feat.qualifiers["gene"][0]
                        elif "locus_tag" in feat.qualifiers and feat.qualifiers["locus_tag"][0] in feature_names:
                            fname = feat.qualifiers["locus_tag"][0]
                        else:
                            fname = None
                    else:
                        if "gene" in feat.qualifiers:
                            fname = feat.qualifiers["gene"][0]
                        else:
                            fname = feat.qualifiers["locus_tag"][0]
                    if fname:
                        features[fname] = feat

            if feature_names is not None:
                for fe in feature_names:
                    if fe not in features:
                        print("Couldn't find gene {} in GFF or GenBank file".format(fe))

    else:
        from Bio import SeqIO
        for feat in SeqIO.read(reference, 'genbank').features:
            if feat.type=='CDS':
                if "locus_tag" in feat.qualifiers:
                    fname = feat.qualifiers["locus_tag"][0]
                    if feature_names is None or fname in feature_names:
                        features[fname] = feat
                elif "gene" in feat.qualifiers:
                    fname = feat.qualifiers["gene"][0]
                    if feature_names is None or fname in feature_names:
                        features[fname] = feat
            elif feat.type=='source': #read 'nuc' as well for annotations - need start/end of whole!
                features['nuc'] = feat

    return features

def read_config(fname):
    if not (fname and os.path.isfile(fname)):
        print("ERROR: config file %s not found."%fname)
        return defaultdict(dict)

    try:
        with open(fname, 'rb') as ifile:
            config = json.load(ifile)
    except json.decoder.JSONDecodeError as err:
        print("FATAL ERROR:")
        print("\tCouldn't parse the JSON file {}".format(fname))
        print("\tError message: '{}'".format(err.msg))
        print("\tLine number: '{}'".format(err.lineno))
        print("\tColumn number: '{}'".format(err.colno))
        print("\tYou must correct this file in order to proceed.")
        sys.exit(2)

    return config

def read_lat_longs(overrides=None, use_defaults=True):
    coordinates = {}
    # TODO: make parsing of tsv files more robust while allow for whitespace delimiting for backwards compatibility
    def add_line_to_coordinates(line):
        if line.startswith('#') or line.strip() == "":
            return
        fields = line.strip().split() if not '\t' in line else line.strip().split('\t')
        if len(fields) == 4:
            geo_field, loc = fields[0].lower(), fields[1].lower()
            lat, long = float(fields[2]), float(fields[3])
            coordinates[(geo_field, loc)] = {
                "latitude": lat,
                "longitude": long
            }
        else:
            print("WARNING: geo-coordinate file contains invalid line. Please make sure not to mix tabs and spaces as delimiters (use only tabs):",line)
    if use_defaults:
        with resource_stream(__package__, "data/lat_longs.tsv") as stream:
            with TextIOWrapper(stream, "utf-8") as defaults:
                for line in defaults:
                    add_line_to_coordinates(line)
    if overrides:
        if os.path.isfile(overrides):
            with open(overrides, encoding='utf-8') as ifile:
                for line in ifile:
                    add_line_to_coordinates(line)
        else:
            print("WARNING: input lat/long file %s not found." % overrides)
    return coordinates

def read_colors(overrides=None, use_defaults=True):
    return ColorParser(mapping_filename=overrides, use_defaults=use_defaults).mapping

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

shquote = shlex.quote

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


def first_line(text):
    """
    Returns the first line of the given text, ignoring leading and trailing
    whitespace.
    """
    return text.strip().splitlines()[0]


def available_cpu_cores(fallback: int = 1) -> int:
    """
    Returns the number (an int) of CPU cores available to this **process**, if
    determinable, otherwise the number of CPU cores available to the
    **computer**, if determinable, otherwise the *fallback* number (which
    defaults to 1).
    """
    try:
        # Note that this is the correct function to use, not os.cpu_count(), as
        # described in the latter's documentation.
        #
        # The reason, which the documentation does not detail, is that
        # processes may be pinned or restricted to certain CPUs by setting
        # their "affinity".  This is not typical except in high-performance
        # computing environments, but if it is done, then a computer with say
        # 24 total cores may only allow our process to use 12.  If we tried to
        # naively use all 24, we'd end up with two threads across the 12 cores.
        # This would degrade performance rather than improve it!
        return len(os.sched_getaffinity(0))
    except:
        # cpu_count() returns None if the value is indeterminable.
        return os.cpu_count() or fallback


def nthreads_value(value):
    """
    Argument value validation and casting function for --nthreads.
    """

    if value.lower() == 'auto':
        return available_cpu_cores()

    try:
        return int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("'%s' is not an integer or the word 'auto'" % value) from None


def get_parent_name_by_child_name_for_tree(tree):
    '''
    Return dictionary mapping child node names to parent node names
    '''
    parents = {}
    for clade in tree.find_clades(order='level'):
        for child in clade:
            parents[child.name] = clade.name
    return parents


def annotate_parents_for_tree(tree):
    """Annotate each node in the given tree with its parent.

    >>> import io
    >>> tree = Bio.Phylo.read(io.StringIO("(A, (B, C))"), "newick")
    >>> not any([hasattr(node, "parent") for node in tree.find_clades()])
    True
    >>> tree = annotate_parents_for_tree(tree)
    >>> tree.root.parent is None
    True
    >>> all([hasattr(node, "parent") for node in tree.find_clades()])
    True
    """
    tree.root.parent = None
    for node in tree.find_clades(order="level"):
        for child in node.clades:
            child.parent = node

    # Return the tree.
    return tree


def json_to_tree(json_dict, root=True):
    """Returns a Bio.Phylo tree corresponding to the given JSON dictionary exported
    by `tree_to_json`.

    Assigns links back to parent nodes for the root of the tree.

    Test opening a JSON from augur export v1.

    >>> import json
    >>> json_fh = open("tests/data/json_tree_to_nexus/flu_h3n2_ha_3y_tree.json", "r")
    >>> json_dict = json.load(json_fh)
    >>> tree = json_to_tree(json_dict)
    >>> tree.name
    'NODE_0002020'
    >>> len(tree.clades)
    2
    >>> tree.clades[0].name
    'NODE_0001489'
    >>> hasattr(tree, "attr")
    True
    >>> "dTiter" in tree.attr
    True
    >>> tree.clades[0].parent.name
    'NODE_0002020'
    >>> tree.clades[0].branch_length > 0
    True

    Test opening a JSON from augur export v2.

    >>> json_fh = open("tests/data/zika.json", "r")
    >>> json_dict = json.load(json_fh)
    >>> tree = json_to_tree(json_dict)
    >>> hasattr(tree, "name")
    True
    >>> len(tree.clades) > 0
    True
    >>> tree.clades[0].branch_length > 0
    True
    """
    # Check for v2 JSON which has combined metadata and tree data.
    if root and "meta" in json_dict and "tree" in json_dict:
        json_dict = json_dict["tree"]

    node = Bio.Phylo.Newick.Clade()

    # v1 and v2 JSONs use different keys for strain names.
    if "name" in json_dict:
        node.name = json_dict["name"]
    else:
        node.name = json_dict["strain"]

    if "children" in json_dict:
        # Recursively add children to the current node.
        node.clades = [json_to_tree(child, root=False) for child in json_dict["children"]]

    # Assign all non-children attributes.
    for attr, value in json_dict.items():
        if attr != "children":
            setattr(node, attr, value)

    # Only v1 JSONs support a single `attr` attribute.
    if hasattr(node, "attr"):
        node.numdate = node.attr.get("num_date")
        node.branch_length = node.attr.get("div")

        if "translations" in node.attr:
            node.translations = node.attr["translations"]
    elif hasattr(node, "node_attrs"):
        node.branch_length = node.node_attrs.get("div")

    if root:
        node = annotate_parents_for_tree(node)

    return node

def get_augur_version():
    """
    Returns a string of the current augur version.
    """
    return __version__


def read_bed_file(bed_file):
    """Read a BED file and return a list of excluded sites.

    Note: This function assumes the given file is a BED file. On parsing
    failures, it will attempt to skip the first line and retry, but no
    other error checking is attempted. Incorrectly formatted files will
    raise errors.

    Parameters
    ----------
    bed_file : str
        Path to the BED file

    Returns:
    --------
    list[int]:
        Sorted list of unique zero-indexed sites
    """
    mask_sites = []
    try:
        bed = pd.read_csv(bed_file, sep='\t', header=None, usecols=[1,2],
                          dtype={1:int,2:int})
    except ValueError:
        # Check if we have a header row. Otherwise, just fail.
        bed = pd.read_csv(bed_file, sep='\t', header=None, usecols=[1,2],
                          dtype={1:int,2:int}, skiprows=1)
        print("Skipped row 1 of %s, assuming it is a header." % bed_file)
    for _, row in bed.iterrows():
        mask_sites.extend(range(row[1], row[2]))
    return sorted(set(mask_sites))

def read_mask_file(mask_file):
    """Read a masking file and return a list of excluded sites.

    Masking files have a single masking site per line, either alone
    or as the second column of a tab-separated file. These sites
    are assumed to be one-indexed, NOT zero-indexed. Incorrectly
    formatted lines will be skipped.

    Parameters
    ----------
    mask_file : str
        Path to the masking file

    Returns:
    --------
    list[int]:
        Sorted list of unique zero-indexed sites
    """
    mask_sites = []
    with open(mask_file, encoding='utf-8') as mf:
        for idx, line in enumerate(l.strip() for l in mf.readlines()):
            if "\t" in line:
                line = line.split("\t")[1]
            try:
                mask_sites.append(int(line) - 1)
            except ValueError as err:
                print("Could not read line %s of %s: '%s' - %s" %
                      (idx, mask_file, line, err), file=sys.stderr)
                raise
    return sorted(set(mask_sites))

def load_mask_sites(mask_file):
    """Load masking sites from either a BED file or a masking file.

    Parameters
    ----------
    mask_file: str
        Path to the BED or masking file

    Returns
    -------
    list[int]
        Sorted list of unique zero-indexed sites
    """
    if mask_file.lower().endswith(".bed"):
        mask_sites = read_bed_file(mask_file)
    else:
        mask_sites = read_mask_file(mask_file)
    print("%d masking sites read from %s" % (len(mask_sites), mask_file))
    return mask_sites

VALID_NUCLEOTIDES = { # http://reverse-complement.com/ambiguity.html
    "A", "G", "C", "T", "U", "N", "R", "Y", "S", "W", "K", "M", "B", "V", "D", "H", "-",
    "a", "g", "c", "t", "u", "n", "r", "y", "s", "w", "k", "m", "b", "v", "d", "h", "-"
}


def read_strains(*files, comment_char="#"):
    """Reads strain names from one or more plain text files and returns the
    set of distinct strains.

    Strain names can be commented with full-line or inline comments. For
    example, the following is a valid strain names file:

        # this is a comment at the top of the file
        strain1  # exclude strain1 because it isn't sequenced properly
        strain2
          # this is an empty line that will be ignored.

    Parameters
    ----------
    files : one or more str
        one or more names of text files with one strain name per line

    Returns
    -------
    set :
        strain names from the given input files

    """
    strains = set()
    for input_file in files:
        with open_file(input_file, 'r') as ifile:
            for line in ifile:
                # Allow comments anywhere in a given line.
                strain_name = line.split(comment_char)[0].strip()
                if len(strain_name) > 0:
                    strains.add(strain_name)

    return strains
