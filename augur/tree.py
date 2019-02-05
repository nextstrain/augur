"""
Build a tree using a variety of methods.
"""

import os
import shutil
import sys
import time
import uuid
import Bio
from Bio import Phylo
import numpy as np
from treetime.vcf_utils import read_vcf
from pathlib import Path

from .utils import run_shell_command, nthreads_value

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
        print("Hint: You can install the missing program using conda or homebrew or apt-get.")
        raise Exception

    return exe


def build_raxml(aln_file, out_file, clean_up=True, nthreads=1):
    '''
    build tree using RAxML with parameters '-f d -m GTRCAT -c 25 -p 235813 -n tre"
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

    call = [raxml,"-T",str(nthreads)," -f d -m GTRCAT -c 25 -p 235813 -n %s -s"%(random_string), aln_file, "> RAxML_log.%s"%(random_string)]
    cmd = " ".join(call)
    print("Building a tree via:\n\t" + cmd +
          "\n\tStamatakis, A: RAxML Version 8: A tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies."
          "\n\tIn Bioinformatics, 2014\n")
    try:
        run_shell_command(cmd, raise_errors = True)
        shutil.copy("RAxML_bestTree.%s"%(random_string), out_file)
        T = Phylo.read(out_file, 'newick')
        if clean_up:
            os.remove("RAxML_bestTree.%s"%(random_string))
            os.remove("RAxML_info.%s"%(random_string))
            os.remove("RAxML_log.%s"%(random_string))
            os.remove("RAxML_parsimonyTree.%s"%(random_string))
            os.remove("RAxML_result.%s"%(random_string))

    except:
        print("TREE BUILDING FAILED, please inspect the raxml.log file\n")
        T=None

    return T

def build_fasttree(aln_file, out_file, clean_up=True, nthreads=1):
    '''
    build tree using fasttree with parameters "-nt"
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

    call = [fasttree, "-nosupport", "-nt", aln_file, "1>", out_file, "2>", log_file]
    cmd = " ".join(call)
    print("Building a tree via:\n\t" + cmd +
          "\n\tPrice et al: FastTree 2 - Approximately Maximum-Likelihood Trees for Large Alignments." +
          "\n\tPLoS ONE 5(3): e9490. https://doi.org/10.1371/journal.pone.0009490\n")
    try:
        run_shell_command(cmd, raise_errors = True, extra_env = extra_env)
        T = Phylo.read(out_file, 'newick')
        if clean_up:
            os.remove(log_file)
    except:
        print("TREE BUILDING FAILED")
        T=None

    return T


def build_iqtree(aln_file, out_file, substitution_model="GTR", clean_up=True, nthreads=1):
    '''
    build tree using IQ-Tree with parameters "-fast"
    arguments:
        aln_file    file name of input aligment
        out_file    file name to write tree to
    '''
    with open(aln_file) as ifile:
        tmp_seqs = ifile.readlines()

    # IQ-tree messes with taxon names. Hence remove offending characters, reinstaniate later
    tmp_aln_file = aln_file.replace(".fasta", "-delim.fasta")
    log_file = tmp_aln_file.replace(".fasta", ".iqtree.log")
    with open(tmp_aln_file, 'w') as ofile:
        for line in tmp_seqs:
            ofile.write(line.replace('/', '_X_X_').replace('|','_Y_Y_'))

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
    fast_opts = [
        "-ninit", "2",
        "-n",     "2",
        "-me",    "0.05"
    ]

    if substitution_model.lower() != "none":
        call = ["iqtree", *fast_opts, "-nt", str(nthreads), "-s", tmp_aln_file,
                "-m", substitution_model, ">", log_file]
    else:
        call = ["iqtree", *fast_opts, "-nt", str(nthreads), "-s", tmp_aln_file, ">", log_file]

    cmd = " ".join(call)

    print("Building a tree via:\n\t" + cmd +
          "\n\tNguyen et al: IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies."
          "\n\tMol. Biol. Evol., 32:268-274. https://doi.org/10.1093/molbev/msu300\n")
    if substitution_model.lower() == "none":
        print("Conducting a model test... see iqtree.log for the result. You can specify this with --substitution-model in future runs.")

    try:
        run_shell_command(cmd, raise_errors = True)
        T = Phylo.read(tmp_aln_file+".treefile", 'newick')
        shutil.copyfile(tmp_aln_file+".treefile", out_file)
        for n in T.get_terminals():
            n.name = n.name.replace('_X_X_','/').replace('_Y_Y_','|')
        #this allows the user to check intermediate output, as tree.nwk will be
        if clean_up:
            #allow user to see chosen model if modeltest was run
            if substitution_model.lower() == 'none':
                shutil.copyfile(log_file, out_file.replace(out_file.split('/')[-1],"iqtree.log"))

            for f in [log_file, tmp_aln_file]:
                if os.path.isfile(f):
                    os.remove(f)

            for ext in [".bionj",".ckp.gz",".iqtree",".log",".mldist",".model.gz",".treefile",".uniqueseq.phy",".model"]:
                if os.path.isfile(tmp_aln_file + ext):
                    os.remove(tmp_aln_file + ext)
    except:
        print("TREE BUILDING FAILED")
        T=None
    return T


def load_excluded_sites(excluded_sites_file):
    """Returns an array of zero-based sites to exclude from a FASTA prior to tree building.

    Parameters
    ----------
    excluded_sites_file : str
        a path to a BED file (with a .bed extension), a tab-delimited DRM file, or a plain text file with one position per line

    Returns
    -------
    ndarray :
        a unique array of positions loaded from the given file
    """
    strip_pos = []
    is_bed_format = False

    if excluded_sites_file is not None:
        # Check for BED file extension.
        if excluded_sites_file.lower().endswith('.bed'):
            import pandas as pd
            is_bed_format = True
            bed = pd.read_csv(excluded_sites_file, sep='\t')
            for index, row in bed.iterrows():
                strip_pos.extend(list(range(row[1], row[2]+1)))
        else:
            # Next, check for DRM-file format or site-per-line format.
            with open(excluded_sites_file, 'r') as ifile:
                line1 = ifile.readline()
                # If the file is tab-delimited, assume it is in DRM-file format.
                if '\t' in line1:
                    strip_pos = [int(line.strip().split('\t')[1]) for line in ifile]
                else:
                    # Finally, fall back to site-per-line format.
                    strip_pos = [int(line.strip()) for line in ifile]
                    if line1.strip():
                        # Add the first line back to the list
                        strip_pos.append(int(line1.strip()))

    strip_pos = np.unique(strip_pos)

    # If the given sites are not in BED format, they are one-based positions and
    # need to be internally adjusted to zero-based positions.
    if not is_bed_format:
        strip_pos = strip_pos - 1

    return strip_pos


def write_out_informative_fasta(compress_seq, alignment, stripFile=None):
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    sequences = compress_seq['sequences']
    ref = compress_seq['reference']
    positions = compress_seq['positions']

    #If want to exclude sites from initial treebuild, read in here
    strip_pos = load_excluded_sites(stripFile)

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
    align = np.rot90(sites)
    seqNamesCorr = list(reversed(seqNames))
    toFasta = [ SeqRecord(id=seqNamesCorr[i], seq=Seq("".join(align[i])), description='') for i in range(len(sequences.keys()))]

    fasta_file = '/'.join(alignment.split('/')[:-1]) + '/informative_sites.fasta'

    #now output this as fasta to read into raxml or iqtree
    SeqIO.write(toFasta, fasta_file, 'fasta')

    #If want a position map, print:
    if printPositionMap:
        with open(fasta_file+".positions.txt", 'w') as the_file:
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
    # Load zero-based excluded sites.
    excluded_sites = load_excluded_sites(excluded_sites_file).tolist()

    # Return the original alignment file, if no excluded sites were found.
    if len(excluded_sites) == 0:
        return alignment_file

    # Load alignment as FASTA generator to prevent loading the whole alignment
    # into memory.
    alignment = Bio.SeqIO.parse(alignment_file, "fasta")

    # Write the masked alignment to disk one record at a time.
    alignment_file_path = Path(alignment_file)
    masked_alignment_file = str(alignment_file_path.parent / ("masked_%s" % alignment_file_path.name))
    with open(masked_alignment_file, "w") as oh:
        for record in alignment:
            # Convert to a mutable sequence to enable masking with Ns.
            sequence = record.seq.tomutable()

            # Replace all excluded sites with Ns.
            for site in excluded_sites:
                sequence[site] = "N"

            record.seq = sequence
            Bio.SeqIO.write(record, oh, "fasta")

    # Return the new alignment FASTA filename.
    return masked_alignment_file


def register_arguments(parser):
    parser.add_argument('--alignment', '-a', required=True, help="alignment in fasta or VCF format")
    parser.add_argument('--method', default='iqtree', choices=["fasttree", "raxml", "iqtree"], help="tree builder to use")
    parser.add_argument('--output', '-o', type=str, help='file name to write tree to')
    parser.add_argument('--substitution-model', default="GTR", choices=["HKY", "GTR", "HKY+G", "GTR+G"],
                                help='substitution model to use. Specify \'none\' to run ModelTest. Currently, only available for IQTREE.')
    parser.add_argument('--nthreads', type=nthreads_value, default=1,
                                help="number of threads to use; specifying the value 'auto' will cause the number of available CPU cores on your system, if determinable, to be used")
    parser.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the VCF was mapped to')
    parser.add_argument('--exclude-sites', type=str, help='file name of one-based sites to exclude for raw tree building (BED format in .bed files, DRM format in tab-delimited files, or one position per line)')


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

    if args.substitution_model and not args.method=='iqtree':
        print("Cannot specify model unless using IQTree. Model specification ignored.")

    if args.method=='raxml':
        T = build_raxml(fasta, tree_fname, nthreads=args.nthreads)
    elif args.method=='iqtree':
        T = build_iqtree(fasta, tree_fname, args.substitution_model, nthreads=args.nthreads)
    elif args.method=='fasttree':
        T = build_fasttree(fasta, tree_fname, nthreads=args.nthreads)
    else:
        print("ERROR: unknown tree builder provided to --method: %s" % args.method, file = sys.stderr)
        return 1

    end = time.time()
    print("Building original tree took {} seconds".format(str(end-start)))

    if T:
        import json
        tree_success = Phylo.write(T, tree_fname, 'newick', format_branch_length='%1.8f')
    else:
        return 1
