import argparse
import Bio
import Bio.Phylo
import os, json, sys
import pandas as pd
import subprocess
import shlex
from treetime.utils import numeric_date
from collections import defaultdict
from pkg_resources import resource_stream
from io import TextIOWrapper
from textwrap import dedent
from .__version__ import __version__
import packaging.version as packaging_version
from .validate import validate, ValidateError, load_json_schema

class AugurException(Exception):
    pass

def myopen(fname, mode):
    if fname.endswith('.gz'):
        import gzip
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)

def get_json_name(args, default=None):
    if args.output:
        print("WARNING: the --output flag will be deprecated in the next major augur release. Use --output-node-data instead.", file=sys.stderr)
        return args.output
    elif args.output_node_data:
        return args.output_node_data
    else:
        if default:
            print("WARNING: no name for the output file was specified. Writing results to %s."%default, file=sys.stderr)
            return default
        else:
            raise ValueError("Please specify a name for the JSON file containing the results.")


def ambiguous_date_to_date_range(mydate, fmt, min_max_year=None):
    from datetime import datetime
    sep = fmt.split('%')[1][-1]
    min_date, max_date = {}, {}
    today = datetime.today().date()

    for val, field  in zip(mydate.split(sep), fmt.split(sep+'%')):
        f = 'year' if 'y' in field.lower() else ('day' if 'd' in field.lower() else 'month')
        if 'XX' in val:
            if f=='year':
                if min_max_year:
                    min_date[f]=min_max_year[0]
                    if len(min_max_year)>1:
                        max_date[f]=min_max_year[1]
                    elif len(min_max_year)==1:
                        max_date[f]=4000 #will be replaced by 'today' below.
                else:
                    return None, None
            elif f=='month':
                min_date[f]=1
                max_date[f]=12
            elif f=='day':
                min_date[f]=1
                max_date[f]=31
        else:
            min_date[f]=int(val)
            max_date[f]=int(val)
    max_date['day'] = min(max_date['day'], 31 if max_date['month'] in [1,3,5,7,8,10,12]
                                           else 28 if max_date['month']==2 else 30)
    lower_bound = datetime(year=min_date['year'], month=min_date['month'], day=min_date['day']).date()
    upper_bound = datetime(year=max_date['year'], month=max_date['month'], day=max_date['day']).date()
    return (lower_bound, upper_bound if upper_bound<today else today)

def read_metadata(fname):
    if not fname:
        print("ERROR: read_metadata called without a filename")
        return {}, []
    if os.path.isfile(fname):
        try:
            metadata = pd.read_csv(fname, sep='\t' if fname[-3:]=='tsv' else ',',
                                    skipinitialspace=True).fillna('')
        except pd.errors.ParserError as e:
            print("Error reading metadata file {}".format(fname))
            print(e)
            sys.exit(2)
        meta_dict = {}
        for ii, val in metadata.iterrows():
            if hasattr(val, "strain"):
                if val.strain in meta_dict:
                    raise ValueError("Duplicate strain '{}'".format(val.strain))
                meta_dict[val.strain] = val.to_dict()
            elif hasattr(val, "name"):
                if val.name in meta_dict:
                    raise ValueError("Duplicate name '{}'".format(val.name))
                meta_dict[val.name] = val.to_dict()
            else:
                print("ERROR: meta data file needs 'name' or 'strain' column")

        return meta_dict, list(metadata.columns)
    else:
        print("ERROR: meta data file ({}) does not exist".format(fname))
        return {}, []


def get_numerical_dates(meta_dict, name_col = None, date_col='date', fmt=None, min_max_year=None):
    if fmt:
        from datetime import datetime
        numerical_dates = {}
        for k,m in meta_dict.items():
            v = m[date_col]
            if type(v)!=str:
                print("WARNING: %s has an invalid data string:"%k,v)
                continue
            elif 'XX' in v:
                ambig_date = ambiguous_date_to_date_range(v, fmt, min_max_year)
                if ambig_date is None or None in ambig_date:
                    numerical_dates[k] = [None, None] #don't send to numeric_date or will be set to today
                else:
                    numerical_dates[k] = [numeric_date(d) for d in ambig_date]
            else:
                try:
                    numerical_dates[k] = numeric_date(datetime.strptime(v, fmt))
                except:
                    numerical_dates[k] = None
    else:
        numerical_dates = {k:float(v) for k,v in meta_dict.items()}

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
    """
    parses one or more "node-data" JSON files and combines them using custom logic.
    Will exit with a (hopefully) helpful message if errors are detected.

    For each JSON, we expect the top-level key "nodes" to be a dict.
    Generated-by fields will not be included in the returned dict of this function.
    """
    if isinstance(fnames, str):
        fnames = [fnames]
    node_data = {"nodes": {}}
    for fname in fnames:
        if os.path.isfile(fname):
            with open(fname) as jfile:
                tmp_data = json.load(jfile)
            if tmp_data.get("annotations"):
                try:
                    validate(tmp_data.get("annotations"), load_json_schema("schema-annotations.json"), fname)
                except ValidateError as err:
                    print("{} contains an `annotations` block of an invalid JSON format. "
                        "Was it produced by different version of augur the one you are currently using ({})? "
                        "Please check the script / program which produced that JSON file.".format(fname, get_augur_version()))
                    print(err)
                    sys.exit(2)
            try:
                for k,v in tmp_data.items():
                    if k=="nodes":
                        if not isinstance(v, dict):
                            raise AugurException("\"nodes\" key in {} is not a dictionary. Please check the formatting of this JSON!".format(fname))
                        for n,nv in v.items():
                            if n in node_data["nodes"]:
                                node_data["nodes"][n].update(nv)
                            else:
                                node_data["nodes"][n] = nv
                    elif k=="generated_by":
                        # Note that this key is _not_ part of the dict returned from this fn.
                        if v.get("program") == "augur" and not is_augur_version_compatable(v.get("version")):
                            # check that the augur version, if provided, is compatible.
                            # ignore version checking of non-augur produced JSONs
                            raise AugurException("Augur version incompatability detected -- the JSON {} was generated by augur version {} which is "
                                "incompatable with the current augur version ({}). We suggest you rerun the pipeline using the current "
                                "version of augur.".format(fname, v.get("version"), get_augur_version()))
                    elif k in node_data:
                        # Behavior as of 2019-11-07 is to do a top-level merge
                        # of dictionaries. If the value is not a dictionary, we
                        # now have a fatal error with a nice message (note that
                        # before 2019-11-07 this was an unhandled error).
                        # This should be revisited in the future. TODO.
                        if isinstance(node_data[k], dict) and isinstance(v, dict):
                            node_data[k].update(v)
                        else:
                            raise AugurException("\"{}\" key found in multiple JSONs. This is not currently handled by augur, "
                                "unless all values are dictionaries. "
                                "Please check the source of these JSONs.".format(k))
                    else:
                        node_data[k]=v
            except AugurException as e:
                print(e)
                sys.exit(2)
        else:
            print("ERROR: node_data JSON file %s not found. Attempting to proceed without it."%fname)

    if tree and os.path.isfile(tree):
        try:
            T = Bio.Phylo.read(tree, 'newick')
        except:
            print("Failed to read tree from file "+tree, file=sys.stderr)
        else:
            tree_node_names = set([l.name for l in T.find_clades()])
            meta_node_names = set(node_data["nodes"].keys())
            if tree_node_names!=meta_node_names:
                print("Names of nodes (including internal nodes) of tree %s don't"
                    " match node names in the node data files."%tree, file=sys.stderr)
    return node_data


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

    with open(file_name, 'w') as handle:
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

        with open(reference) as in_handle:
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
        with open(fname) as ifile:
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
            with open(overrides) as ifile:
                for line in ifile:
                    add_line_to_coordinates(line)
        else:
            print("WARNING: input lat/long file %s not found." % overrides)
    return coordinates

def read_colors(overrides=None, use_defaults=True):
    colors = {}
    # TODO: make parsing of tsv files more robust while allow for whitespace delimiting for backwards compatibility
    def add_line(line):
        if line.startswith('#'):
            return
        fields = line.strip().split() if not '\t' in line else line.strip().split('\t')
        if not fields:
            return # blank lines
        if len(fields) != 3:
            print("WARNING: Color map file contains invalid line. Please make sure not to mix tabs and spaces as delimiters (use only tabs):",line)
            return
        trait, trait_value, hex_code = fields[0].lower(), fields[1].lower(), fields[2]
        if not hex_code.startswith("#") or len(hex_code) != 7:
            print("WARNING: Color map file contained this invalid hex code: ", hex_code)
            return
        # If was already added, delete entirely so order can change to order in user-specified file
        # (even though dicts shouldn't be relied on to have order)
        if (trait, trait_value) in colors:
            del colors[(trait, trait_value)]
        colors[(trait, trait_value)] = hex_code


    if use_defaults:
        with resource_stream(__package__, "data/colors.tsv") as stream:
            with TextIOWrapper(stream, "utf-8") as defaults:
                for line in defaults:
                    add_line(line)

    if overrides:
        if os.path.isfile(overrides):
            with open(overrides) as fh:
                for line in fh:
                    add_line(line)
        else:
            print("WARNING: Couldn't open color definitions file {}.".format(overrides))
    color_map = defaultdict(list)
    for (trait, trait_value), hex_code in colors.items():
        color_map[trait].append((trait_value, hex_code))

    return color_map

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
    with open(vcf_file_name, 'w') as the_file:
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
    with open(ref_file_name, 'w') as the_file:
        the_file.write("\n".join(refWrite))

    with open(vcf_file_name, 'a') as the_file:
        the_file.write("\n".join(vcfWrite))

    if vcf_file_name.lower().endswith('.gz'):
        import os
        #must temporarily remove .gz ending, or gzip won't zip it!
        os.rename(vcf_file_name, vcf_file_name[:-3])
        call = ["gzip", vcf_file_name[:-3]]
        run_shell_command(" ".join(call), raise_errors = True)

shquote = shlex.quote

def run_shell_command(cmd, raise_errors = False, extra_env = None):
    """
    Run the given command string via Bash with error checking.

    Returns True if the command exits normally.  Returns False if the command
    exits with failure and "raise_errors" is False (the default).  When
    "raise_errors" is True, exceptions are rethrown.

    If an *extra_env* mapping is passed, the provided keys and values are
    overlayed onto the default subprocess environment.
    """
    env = os.environ.copy()

    if extra_env:
        env.update(extra_env)

    shargs = ['-c', "set -euo pipefail; " + cmd]

    if os.name == 'posix':
        shellexec = ['/bin/bash']
    else:
        # We try best effort on other systems. For now that means nt/java.
        shellexec = ['env', 'bash']

    try:
        # Use check_call() instead of run() since the latter was added only in Python 3.5.
        subprocess.check_output(
            shellexec + shargs,
            shell = False,
            stderr = subprocess.STDOUT,
            env = env)

    except subprocess.CalledProcessError as error:
        print_error(
            "{out}\nshell exited {rc} when running: {cmd}{extra}",
            out = error.output,
            rc  = error.returncode,
            cmd = cmd,
            extra = "\nAre you sure this program is installed?" if error.returncode==127 else "",
        )
        if raise_errors:
            raise
        else:
            return False

    except FileNotFoundError as error:
        print_error(
            """
            Unable to run shell commands using {shell}!

            Augur requires {shell} to be installed.  Please open an issue on GitHub
            <https://github.com/nextstrain/augur/issues/new> if you need assistance.
            """,
            shell = ' and '.join(shellexec)
        )
        if raise_errors:
            raise
        else:
            return False

    else:
        return True


def print_error(message, **kwargs):
    """
    Formats *message* with *kwargs* using :meth:`str.format` and
    :func:`textwrap.dedent` and uses it to print an error message to
    ``sys.stderr``.
    """
    print("\nERROR: " + dedent(message.format(**kwargs)).lstrip("\n")+"\n", file = sys.stderr)


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

def is_augur_version_compatable(version):
    """
    Checks if the provided **version** is the same major version
    as the currently running version of augur.

    Parameters
    ----------
    version : str
        version to check against the current version

    Returns
    -------
    Bool

    """
    current_version = packaging_version.parse(get_augur_version())
    this_version = packaging_version.parse(version)
    return this_version.release[0] == current_version.release[0]
