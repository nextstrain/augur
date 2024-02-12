"""
Build a tree using a variety of methods.
"""

import os
import shlex
import shutil
import sys
import time
import uuid
import Bio
from Bio.Seq import MutableSeq
from Bio import Phylo
import numpy as np
from treetime.vcf_utils import read_vcf
from pathlib import Path

from .errors import AugurError
from .io.file import open_file
from .io.sequences import read_sequences
from .io.shell_command_runner import run_shell_command
from .io.vcf import shquote
from .utils import nthreads_value, load_mask_sites

DEFAULT_ARGS = {
    "fasttree": "-nt -nosupport",
    "raxml": "-f d -m GTRCAT -c 25 -p 235813",
    # For compat with older versions of iqtree, we avoid the newish -fast
    # option alias and instead spell out its component parts:
    #
    #     -ninit 2
    #     -n 2
    #     -me 0.05
    #
    # This may need to be updated in the future if we want to stay in lock-step
    # with -fast, although there's probably no particular reason we have to.
    # Refer to the handling of -fast in utils/tools.cpp:
    #   https://github.com/Cibiv/IQ-TREE/blob/44753aba/utils/tools.cpp#L2926-L2936
    # Increasing threads (nt) can cause IQtree to run longer, hence use AUTO by default
    # Auto makes IQtree chose the optimal number of threads
    # Redo prevents IQtree errors when a run was aborted and restarted
    "iqtree": "-ninit 2 -n 2 -me 0.05 -nt AUTO -redo",
}

# IQ-TREE only; see usage below
DEFAULT_SUBSTITUTION_MODEL = "GTR"

class ConflictingArgumentsException(Exception):
    """Exception when user-provided tree builder arguments conflict with the
    requested tree builder's hardcoded defaults (e.g., the path to the
    alignment, etc.).

    """
    pass

def check_conflicting_args(tree_builder_args, defaults):
    """Checks the given user-provided tree builder arguments for hardcoded default
    arguments and raise an exception with a list of any that are found.

    Arguments
    ---------
    tree_builder_args : str
        User-provided tree builder arguments
    defaults : list or tuple
        List of hardcoded default arguments (e.g., ['-nt'])

    Raises
    ------
    ConflictingArgumentsException
        When any user-provided arguments match those in the defaults.

    Examples
    --------
    >>> defaults = ("-ntmax", "-m", "-s")
    >>> check_conflicting_args("-czb -n 2", defaults)
    >>> check_conflicting_args("-czb -ntmax 2", defaults)
    Traceback (most recent call last):
        ...
    augur.tree.ConflictingArgumentsException: The following tree builder arguments conflict with hardcoded defaults. Remove these arguments and try again: -ntmax

    """
    # Parse tree builder argument string into a list of shell arguments. This
    # allows us to search a list of arguments instead of a string where we might
    # find partial prefix matches to some of the given defaults.
    tree_builder_args_list = shlex.split(tree_builder_args)
    conflicting_args = [
        default
        for default in defaults
        if default in tree_builder_args_list
    ]

    if len(conflicting_args) > 0:
        raise ConflictingArgumentsException(
            f"The following tree builder arguments conflict with hardcoded defaults. Remove these arguments and try again: {', '.join(conflicting_args)}"
        )

def find_executable(names, default = None):
    """
    Return the path to the first executable found in PATH from the given list
    of names.

    Raises a (hopefully helpful) error if no executable is found.  Provide a
    value for the "default" parameter to instead return a value.
    """
    exe = next(filter(shutil.which, names), default)

    if exe is None:
        print("Unable to find any of %s in PATH=%s" % (names, os.environ["PATH"]))
        print("\nHint: You can install the missing program using conda or homebrew or apt-get.\n")
        raise Exception

    return exe


def build_raxml(aln_file, out_file, clean_up=True, nthreads=1, tree_builder_args=None):
    '''
    build tree using RAxML
    '''

    raxml = find_executable([
        # Users who symlink/install as "raxml" can pick a specific version,
        # otherwise we have our own search order based on expected parallelism.
        # The variants we look for are not exhaustive, but based on what's
        # provided by conda and Ubuntu's raxml packages.  This may want to be
        # adjusted in the future depending on use.
        "raxml",
        "raxmlHPC-PTHREADS-AVX2",
        "raxmlHPC-PTHREADS-AVX",
        "raxmlHPC-PTHREADS-SSE3",
        "raxmlHPC-PTHREADS",
        "raxmlHPC-AVX2",
        "raxmlHPC-AVX",
        "raxmlHPC-SSE3",
        "raxmlHPC",
    ])

    # RAxML outputs files appended with this random string:
    # RAxML_bestTree.4ed91a, RAxML_info.4ed91a, RAxML_parsimonyTree.4ed91a, RAxML_result.4ed91a
    random_string = uuid.uuid4().hex[0:6]

    # Check tree builder arguments for conflicts with hardcoded defaults.
    check_conflicting_args(tree_builder_args, ("-T", "-n", "-s"))

    call = [raxml,"-T",str(nthreads), "-n %s -s"%(random_string), shquote(aln_file), tree_builder_args, "> RAxML_log.%s"%(random_string)]
    cmd = " ".join(call)
    print("Building a tree via:\n\t" + cmd +
          "\n\tStamatakis, A: RAxML Version 8: A tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies."
          "\n\tIn Bioinformatics, 2014\n")
    try:
        run_shell_command(cmd, raise_errors = True)

        # When bootstrapping is enabled by custom tree_builder_args, RAxML
        # outputs the tree with support values to a different file,
        # "bipartitions".  If it exists, then use it; otherwise use "bestTree".
        possible_tree_files = [
            Path(f"RAxML_bipartitions.{random_string}"),    # best-scoring ML tree with support values
            Path(f"RAxML_bestTree.{random_string}"),        # best-scoring ML tree
        ]

        tree = next((t for t in possible_tree_files if t.exists()), None)

        if not tree:
            raise AugurError(f"No RAxML output tree files found; looked for: {', '.join(map(repr, map(str, possible_tree_files)))}")

        shutil.copy(str(tree), out_file)
        T = Phylo.read(out_file, 'newick')

        if clean_up:
            for f in Path().glob(f"RAxML_*.{random_string}"):
                f.unlink()

    except Exception as error:
        print("ERROR: TREE BUILDING FAILED")
        print(f"ERROR: {error}")
        if os.path.isfile("RAxML_log.%s"%(random_string)):
            print("Please see the log file for more details: {}".format("RAxML_log.%s"%(random_string)))
        T=None

    return T


def build_fasttree(aln_file, out_file, clean_up=True, nthreads=1, tree_builder_args=None):
    '''
    build tree using fasttree
    '''
    log_file = out_file + ".log"

    fasttree = find_executable([
        # Search order is based on expected parallelism and accuracy
        # (double-precision versions).
        "FastTreeDblMP",
        "FastTreeDbl",
        "FastTreeMP",
        "fasttreeMP",   # Ubuntu lowercases
        "FastTree",
        "fasttree"
    ])

    # By default FastTree with OpenMP support uses all available cores.
    # However, it respects the standard OpenMP environment variable controlling
    # this as described at <http://www.microbesonline.org/fasttree/#OpenMP>.
    #
    # We always set it, regardless of it the found FastTree executable contains
    # "MP" in the name, because the generic "FastTree" or "fasttree" variants
    # might be OpenMP-enabled too.
    extra_env = {
        "OMP_NUM_THREADS": str(nthreads),
    }

    call = [fasttree, tree_builder_args, shquote(aln_file), "1>", shquote(out_file), "2>", shquote(log_file)]
    cmd = " ".join(call)
    print("Building a tree via:\n\t" + cmd +
          "\n\tPrice et al: FastTree 2 - Approximately Maximum-Likelihood Trees for Large Alignments." +
          "\n\tPLoS ONE 5(3): e9490. https://doi.org/10.1371/journal.pone.0009490\n")
    try:
        run_shell_command(cmd, raise_errors = True, extra_env = extra_env)
        T = Phylo.read(out_file, 'newick')
    except Exception as error:
        print("ERROR: TREE BUILDING FAILED")
        print(f"ERROR: {error}")
        if os.path.isfile(log_file):
            print("Please see the log file for more details: {}".format(log_file))
        T=None

    return T


def build_iqtree(aln_file, out_file, substitution_model="GTR", clean_up=True, nthreads=1, tree_builder_args=None):
    '''
    build tree using IQ-Tree
    arguments:
        aln_file    file name of input aligment
        out_file    file name to write tree to
    '''
    iqtree = find_executable([
        "iqtree2",
        "iqtree"
    ])
    # create a dictionary for characters that IQ-tree changes.
    # we remove those prior to tree-building and reinstantiate later
    def random_string(n):
        from string import ascii_uppercase as letters
        return "".join([letters[i] for i in np.random.randint(len(letters), size=n)])
    prefix = "DELIM"
    escape_dict = {c:f'_{prefix}-{random_string(20)}_' for c in '/|()*'}
    reverse_escape_dict = {v:k for k,v in escape_dict.items()}


    # IQ-tree messes with taxon names. Hence remove offending characters, reinstaniate later
    tmp_aln_file = str(Path(aln_file).with_name(Path(aln_file).stem + "-delim.fasta"))
    log_file = str(Path(tmp_aln_file).with_suffix(".iqtree.log"))
    num_seqs = 0
    with open_file(tmp_aln_file, 'w') as ofile, open_file(aln_file) as ifile:
        for line in ifile:
            tmp_line = line
            if line.startswith(">"):
                num_seqs += 1
                for c,v in escape_dict.items():
                    tmp_line = tmp_line.replace(c,v)

            ofile.write(tmp_line)

    # Check tree builder arguments for conflicts with hardcoded defaults.
    check_conflicting_args(tree_builder_args, ("-ntmax", "-s", "-m"))

    if substitution_model.lower() != "auto":
        call = [iqtree, "-ntmax", str(nthreads), "-s", shquote(tmp_aln_file),
                "-m", shquote(substitution_model), tree_builder_args, ">", shquote(log_file)]
    else:
        call = [iqtree, "-ntmax", str(nthreads), "-s", shquote(tmp_aln_file), tree_builder_args, ">", shquote(log_file)]

    cmd = " ".join(call)

    print("Building a tree via:\n\t" + cmd +
          "\n\tNguyen et al: IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies."
          "\n\tMol. Biol. Evol., 32:268-274. https://doi.org/10.1093/molbev/msu300\n")
    if substitution_model.lower() == "auto":
        print(f"Conducting a model test... see '{shquote(log_file)}' for the result. You can specify this with --substitution-model in future runs.")

    try:
        run_shell_command(cmd, raise_errors = True)
        T = Phylo.read(tmp_aln_file+".treefile", 'newick')
        shutil.copyfile(tmp_aln_file+".treefile", out_file)
        for n in T.find_clades(terminal=True):
            tmp_name = n.name
            for v,c in reverse_escape_dict.items():
                tmp_name = tmp_name.replace(v,c)
            n.name = tmp_name
        #this allows the user to check intermediate output, as tree.nwk will be
        if clean_up:
            if os.path.isfile(tmp_aln_file):
                os.remove(tmp_aln_file)

            for ext in [".bionj",".ckp.gz",".iqtree",".mldist",".model.gz",".treefile",".uniqueseq.phy",".model"]:
                if os.path.isfile(tmp_aln_file + ext):
                    os.remove(tmp_aln_file + ext)
    except Exception as error:
        print("ERROR: TREE BUILDING FAILED")
        print(f"ERROR: {error}")
        if os.path.isfile(log_file):
            print("Please see the log file for more details: {}".format(log_file))
        T=None
    return T


def write_out_informative_fasta(compress_seq, alignment, stripFile=None):
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    sequences = compress_seq['sequences']
    ref = compress_seq['reference']
    positions = compress_seq['positions']

    #If want to exclude sites from initial treebuild, read in here
    strip_pos = load_mask_sites(stripFile) if stripFile else []

    #Get sequence names
    seqNames = list(sequences.keys())

    #Check non-ref sites to see if informative
    printPositionMap = False    #If true, prints file mapping Fasta position to real position
    sites = []
    pos = []

    for key in positions:
        if key not in strip_pos:
            pattern = []
            for k in sequences.keys():
                #looping try/except is faster than list comprehension
                try:
                    pattern.append(sequences[k][key])
                except KeyError:
                    pattern.append(ref[key])
            origPattern = list(pattern)
            if '-' in pattern or 'N' in pattern:
                #remove gaps/Ns to see if otherwise informative
                pattern = [value for value in origPattern if value != '-' and value != 'N']
            un = np.unique(pattern, return_counts=True)
            #If not all - or N, not all same base, and >1 differing base, append
            if len(un[0])!=0 and len(un[0])!=1 and not (len(un[0])==2 and min(un[1])==1):
                sites.append(origPattern)
                pos.append("\t".join([str(len(pos)+1),str(key)]))

    #Rotate and convert to SeqRecord
    sites = np.asarray(sites)
    if len(sites.shape)!=2:
        print("ERROR: NO VALID SITES REMAIN AFTER IGNORING UNCALLED SITES")
        raise Exception

    align = np.rot90(sites)
    seqNamesCorr = list(reversed(seqNames))
    toFasta = [ SeqRecord(id=seqNamesCorr[i], seq=Seq("".join(align[i])), description='') for i in range(len(sequences.keys()))]

    fasta_file = os.path.join(os.path.dirname(alignment), 'informative_sites.fasta')

    #now output this as fasta to read into raxml or iqtree
    SeqIO.write(toFasta, fasta_file, 'fasta')

    #If want a position map, print:
    if printPositionMap:
        with open_file(fasta_file+".positions.txt", 'w') as the_file:
            the_file.write("\n".join(pos))

    return fasta_file


def mask_sites_in_multiple_sequence_alignment(alignment_file, excluded_sites_file):
    """Creates a new multiple sequence alignment FASTA file from which the given
    excluded sites have been removed and returns the filename of the new
    alignment.

    Parameters
    ----------
    alignment_file : str
        path to the original multiple sequence alignment file

    excluded_sites_file : str
        path to a text file containing each nucleotide position to exclude with one position per line

    Returns
    -------
    str
        path to the new FASTA file from which sites have been excluded
    """
    # Read 1-based excluded sites and store as 0-based sites.
    excluded_sites = load_mask_sites(excluded_sites_file)

    # Return the original alignment file, if no excluded sites were found.
    if len(excluded_sites) == 0:
        return alignment_file

    # Load alignment as FASTA generator to prevent loading the whole alignment
    # into memory.
    alignment = read_sequences(alignment_file)

    # Write the masked alignment to disk one record at a time.
    alignment_file_path = Path(alignment_file)
    masked_alignment_file = str(alignment_file_path.parent / ("masked_%s" % alignment_file_path.name))
    with open_file(masked_alignment_file, "w") as oh:
        for record in alignment:
            # Convert to a mutable sequence to enable masking with Ns.
            sequence = MutableSeq(str(record.seq))

            # Replace all excluded sites with Ns.
            for site in excluded_sites:
                sequence[site] = "N"

            record.seq = sequence
            Bio.SeqIO.write(record, oh, "fasta")

    # Return the new alignment FASTA filename.
    return masked_alignment_file


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("tree", help=__doc__)
    parser.add_argument('--alignment', '-a', required=True, help="alignment in fasta or VCF format")
    parser.add_argument('--method', default='iqtree', choices=["fasttree", "raxml", "iqtree"], help="tree builder to use")
    parser.add_argument('--output', '-o', type=str, help='file name to write tree to')
    parser.add_argument('--substitution-model', default=DEFAULT_SUBSTITUTION_MODEL,
                                help='substitution model to use. Specify \'auto\' to run ModelTest. Currently, only available for IQ-TREE.')
    parser.add_argument('--nthreads', type=nthreads_value, default=1,
                                help="maximum number of threads to use; specifying the value 'auto' will cause the number of available CPU cores on your system, if determinable, to be used")
    parser.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the VCF was mapped to')
    parser.add_argument('--exclude-sites', type=str, help='file name of one-based sites to exclude for raw tree building (BED format in .bed files, second column in tab-delimited files, or one position per line)')
    parser.add_argument('--tree-builder-args', type=str, help=f"""arguments to pass to the tree builder either augmenting or overriding the default arguments (except for input alignment path, number of threads, and substitution model).
    Use the assignment operator (e.g., --tree-builder-args="-czb" for IQ-TREE) to avoid unexpected errors.
    FastTree defaults: "{DEFAULT_ARGS['fasttree']}".
    RAxML defaults: "{DEFAULT_ARGS['raxml']}".
    IQ-TREE defaults: "{DEFAULT_ARGS['iqtree']}".
    """)
    parser.add_argument('--override-default-args', action="store_true", help="override default tree builder arguments with the values provided by the user in `--tree-builder-args` instead of augmenting the existing defaults.")

    parser.epilog = """For example, to build a tree with IQ-TREE, use the following format:
    augur tree --method iqtree --alignment <alignment> --substitution-model <model> --output <tree> --tree-builder-args="<extra arguments>"
    """
    return parser

def run(args):
    # check alignment type, set flags, read in if VCF
    is_vcf = False
    ref = None
    if any([args.alignment.lower().endswith(x) for x in ['.vcf', '.vcf.gz']]):
        # Prepare a multiple sequence alignment from the given variants VCF and
        # reference FASTA.
        if not args.vcf_reference:
            print("ERROR: a reference Fasta is required with VCF-format alignments")
            return 1
        compress_seq = read_vcf(args.alignment, args.vcf_reference)
        sequences = compress_seq['sequences']
        ref = compress_seq['reference']
        is_vcf = True
        aln = sequences
    elif args.exclude_sites:
        # Mask excluded sites from the given multiple sequence alignment.
        aln = mask_sites_in_multiple_sequence_alignment(args.alignment, args.exclude_sites)
    else:
        # Use the multiple sequence alignment as is.
        aln = args.alignment

    start = time.time()

    if args.output:
        tree_fname = args.output
    else:
        tree_fname = '.'.join(args.alignment.split('.')[:-1]) + '.nwk'

    # construct reduced alignment if needed
    if is_vcf:
        variable_fasta = write_out_informative_fasta(compress_seq, args.alignment, stripFile=args.exclude_sites)
        fasta = variable_fasta
    else:
        fasta = aln

    if args.method != "iqtree" and args.substitution_model is not DEFAULT_SUBSTITUTION_MODEL:
        print(f"Cannot specify --substitution-model unless using IQ-TREE. Model specification {args.substitution_model!r} ignored.", file=sys.stderr)

    # Allow users to keep default args, override them, or augment them.
    if args.tree_builder_args is None:
        tree_builder_args = DEFAULT_ARGS[args.method]
    elif args.override_default_args:
        tree_builder_args = args.tree_builder_args
    else:
        tree_builder_args = f"{DEFAULT_ARGS[args.method]} {args.tree_builder_args}"

    try:
        if args.method=='raxml':
            T = build_raxml(fasta, tree_fname, nthreads=args.nthreads, tree_builder_args=tree_builder_args)
        elif args.method=='iqtree':
            T = build_iqtree(fasta, tree_fname, args.substitution_model, nthreads=args.nthreads, tree_builder_args=tree_builder_args)
        elif args.method=='fasttree':
            T = build_fasttree(fasta, tree_fname, nthreads=args.nthreads, tree_builder_args=tree_builder_args)
    except ConflictingArgumentsException as error:
        print(f"ERROR:", error, file=sys.stderr)
        return 1

    end = time.time()
    print("\nBuilding original tree took {} seconds".format(str(end-start)))

    if T:
        tree_success = Phylo.write(T, tree_fname, 'newick', format_branch_length='%1.8f')
    else:
        return 1
