"""
Build a tree using a variety of methods.
"""

import os, shutil, sys, time
from Bio import Phylo
from treetime.vcf_utils import read_vcf
import numpy as np
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

    call = [raxml,"-T",str(nthreads)," -f d -m GTRCAT -c 25 -p 235813 -n tre -s",aln_file,"> raxml.log"]
    cmd = " ".join(call)
    print("Building a tree via:\n\t" + cmd +
          "\n\tStamatakis, A: RAxML Version 8: A tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies."
          "\n\tIn Bioinformatics, 2014\n")
    try:
        run_shell_command(cmd, raise_errors = True)
        shutil.copy('RAxML_bestTree.tre', out_file)
        T = Phylo.read(out_file, 'newick')
        if clean_up:
            os.remove('raxml.log')
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
    aln_file = "temp_iqtree.fasta"
    with open(aln_file, 'w') as ofile:
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
        call = ["iqtree", *fast_opts, "-nt", str(nthreads), "-s", aln_file, "-m", substitution_model,
            ">", "iqtree.log"]
    else:
        call = ["iqtree", *fast_opts, "-nt", str(nthreads), "-s", aln_file, ">", "iqtree.log"]

    cmd = " ".join(call)

    print("Building a tree via:\n\t" + cmd +
          "\n\tNguyen et al: IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies."
          "\n\tMol. Biol. Evol., 32:268-274. https://doi.org/10.1093/molbev/msu300\n")
    if substitution_model.lower() == "none":
        print("Conducting a model test... see iqtree.log for the result. You can specify this with --substitution-model in future runs.")

    try:
        run_shell_command(cmd, raise_errors = True)
        T = Phylo.read(aln_file+".treefile", 'newick')
        shutil.copyfile(aln_file+".treefile", out_file)
        #this allows the user to check intermediate output, as tree.nwk will be
        if clean_up:
            #allow user to see chosen model if modeltest was run
            if substitution_model.lower() == 'none':
                shutil.copyfile('iqtree.log', out_file.replace(out_file.split('/')[-1],"iqtree.log"))
            os.remove('iqtree.log')
            os.remove(aln_file)
            for ext in [".bionj",".ckp.gz",".iqtree",".log",".mldist",".model.gz",".treefile",".uniqueseq.phy",".model"]:
                if os.path.isfile(aln_file + ext):
                    os.remove(aln_file + ext)
            for n in T.get_terminals():
                n.name = n.name.replace('_X_X_','/').replace('_Y_Y_','|')
    except:
        print("TREE BUILDING FAILED")
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
    #IF FIND STANDARDIZED DRM FILE FORMAT, IMPLEMENT HERE
    strip_pos = []
    if stripFile:
        if stripFile.lower().endswith('.bed'): #BED format file
            import pandas as pd
            bed = pd.read_csv(stripFile, sep='\t')
            for index, row in bed.iterrows():
                strip_pos.extend(list(range(row[1], row[2]+1)))
        else: #site-per-line format or DRM-file format
            with open(stripFile, 'r') as ifile:
                line1 = ifile.readline()
                if '\t' in line1: #DRM-file format
                    strip_pos = [int(line.strip().split('\t')[1]) for line in ifile]
                else: #site-per-line
                    strip_pos = [int(line.strip()) for line in ifile]
                    strip_pos.append(int(line1.strip())) #add back 1st line
        strip_pos = np.unique(strip_pos)
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


def register_arguments(parser):
    parser.add_argument('--alignment', '-a', required=True, help="alignment in fasta or VCF format")
    parser.add_argument('--method', default='iqtree', choices=["fasttree", "raxml", "iqtree"], help="tree builder to use")
    parser.add_argument('--output', '-o', type=str, help='file name to write tree to')
    parser.add_argument('--substitution-model', default="GTR", choices=["HKY", "GTR", "HKY+G", "GTR+G"],
                                help='substitution model to use. Specify \'none\' to run ModelTest. Currently, only available for IQTREE.')
    parser.add_argument('--nthreads', type=nthreads_value, default=1,
                                help="number of threads to use; specifying the value 'auto' will cause the number of available CPU cores on your system, if determinable, to be used")
    parser.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the VCF was mapped to')
    parser.add_argument('--exclude-sites', type=str, help='file name of sites to exclude for raw tree building (VCF only)')


def run(args):
    # check alignment type, set flags, read in if VCF
    is_vcf = False
    ref = None
    if any([args.alignment.lower().endswith(x) for x in ['.vcf', '.vcf.gz']]):
        if not args.vcf_reference:
            print("ERROR: a reference Fasta is required with VCF-format alignments")
            return 1
        compress_seq = read_vcf(args.alignment, args.vcf_reference)
        sequences = compress_seq['sequences']
        ref = compress_seq['reference']
        is_vcf = True
        aln = sequences
    else:
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
