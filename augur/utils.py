import argparse
import Bio
import Bio.Phylo
import os, json, sys
import pandas as pd
import subprocess
from treetime.utils import numeric_date
from collections import defaultdict
from pkg_resources import resource_stream
from io import TextIOWrapper

def myopen(fname, mode):
    if fname.endswith('.gz'):
        import gzip
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)

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
    if os.path.isfile(fname):
        metadata = pd.read_csv(fname, sep='\t' if fname[-3:]=='tsv' else ',',
                             skipinitialspace=True).fillna('')
        meta_dict = {}
        for ii, val in metadata.iterrows():
            if hasattr(val, "strain"):
                meta_dict[val.strain] = val.to_dict()
            elif hasattr(val, "name"):
                meta_dict[val.name] = val.to_dict()
            else:
                print("ERROR: meta data file needs 'name' or 'strain' column")

        return meta_dict, list(metadata.columns)
    else:
        print("ERROR: file with states does not exist")
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

def read_node_data(fnames, tree=None):
    """parse the "nodes" field of the given JSONs and join the data together"""
    if type(fnames) is str:
        fnames = [fnames]
    node_data = {"nodes": {}}
    for fname in fnames:
        if os.path.isfile(fname):
            with open(fname) as jfile:
                tmp_data = json.load(jfile)
            for k,v in tmp_data.items():
                if k=="nodes":
                    for n,nv in v.items():
                        if n in node_data["nodes"]:
                            node_data["nodes"][n].update(nv)
                        else:
                            node_data["nodes"][n] = nv
                elif k in node_data:
                    # will this recurse into nested dicts?!?!
                    node_data[k].update(v)
                else:
                    node_data[k]=v
        else:
            print("ERROR: node_data JSON file %s not found. Attempting to proceed without it."%fname)

    if tree and os.path.isfile(tree):
        from Bio import Phylo
        try:
            T = Phylo.read(tree, 'newick')
        except:
            print("Failed to read tree from file "+tree, file=sys.stderr)
        else:
            tree_node_names = set([l.name for l in T.find_clades()])
            meta_node_names = set(node_data["nodes"].keys())
            if tree_node_names!=meta_node_names:
                print("Names of nodes (including internal nodes) of tree %s don't"
                      " match node names in the node data files."%tree, file=sys.stderr)

    return node_data


def write_json(data, file_name, indent=1):
    import json, os
    success = False

    #in case auspice folder does not exist yet
    parent_directory = os.path.dirname(file_name)
    if parent_directory and not os.path.exists(parent_directory):
        try:
            os.makedirs(parent_directory)
        except OSError: #Guard against race condition
            if not os.path.isdir(parent_directory):
                raise

    try:
        handle = open(file_name, 'w')
    except IOError:
        raise
    else:
        json.dump(data, handle, indent=indent, sort_keys = True)
        handle.close()
        success=True

    return success


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
    if fname and os.path.isfile(fname):
        with open(fname) as ifile:
            config = json.load(ifile)
    else:
        print("ERROR: config file %s not found."%fname)
        config = defaultdict(dict)

    return config

def read_lat_longs(overrides=None, use_defaults=True):
    coordinates = {}
    def add_line_to_coordinates(line):
        if line.startswith('#'): 
            return
        fields = line.strip().split()
        if len(fields) == 4:
            geo_field, loc = fields[0].lower(), fields[1].lower()
            lat, long = float(fields[2]), float(fields[3])
            coordinates[(geo_field, loc)] = {
                "latitude": lat,
                "longitude": long
            }
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
    def add_line(line):
        if line.startswith('#'):
            return
        fields = line.strip().split()
        if not fields:
            return # blank lines
        if len(fields) != 3:
            print("WARNING: Color map file contained this invalid line: ", line)
            return
        trait, trait_value, hex_code = fields[0].lower(), fields[1].lower(), fields[2]
        if not hex_code.startswith("#") or len(hex_code) != 7:
            print("WARNING: Color map file contained this invalid hex code: ", hex_code)
            return
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


def run_shell_command(cmd, raise_errors = False, extra_env = None):
    """
    Run the given command string via the shell with error checking.

    Returns True if the command exits normally.  Returns False if the command
    exits with failure and "raise_errors" is False (the default).  When
    "raise_errors" is True, exceptions are rethrown.

    If an *extra_env* mapping is passed, the provided keys and values are
    overlayed onto the default subprocess environment.
    """
    env = os.environ.copy()

    if extra_env:
        env.update(extra_env)

    try:
        # Use check_call() instead of run() since the latter was added only in Python 3.5.
        subprocess.check_call(cmd, shell = True, env = env)
    except subprocess.CalledProcessError as error:
        print(
            "ERROR: {program} exited {returncode}, invoked as: {cmd}".format(
                program    = cmd.split()[0],
                returncode = error.returncode,
                cmd        = cmd,
            ),
            file = sys.stderr
        )
        if raise_errors:
            raise
        else:
            return False
    else:
        return True


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
    """
    tree.parent = None
    for node in tree.find_clades(order="level"):
        for child in node.clades:
            child.parent = node

    # Return the tree.
    return tree


def json_to_tree(json_dict, root=True):
    """Returns a Bio.Phylo tree corresponding to the given JSON dictionary exported
    by `tree_to_json`.

    Assigns links back to parent nodes for the root of the tree.

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
    """
    node = Bio.Phylo.Newick.Clade()
    node.name = json_dict["strain"]

    if "children" in json_dict:
        # Recursively add children to the current node.
        node.clades = [json_to_tree(child, root=False) for child in json_dict["children"]]

    # Assign all non-children attributes.
    for attr, value in json_dict.items():
        if attr != "children":
            setattr(node, attr, value)

    node.numdate = node.attr.get("num_date")
    node.branch_length = node.attr.get("div")

    if "translations" in node.attr:
        node.translations = node.attr["translations"]

    if root:
        node = annotate_parents_for_tree(node)

    return node
