import os, shutil, time
from Bio import Phylo
from .utils import read_metadata, get_numerical_dates, write_json
from treetime.vcf_utils import read_vcf, write_vcf
import numpy as np

def build_raxml(aln_file, out_file, clean_up=True, nthreads=2):
    '''
    build tree using RAxML with parameters '-f d -m GTRCAT -c 25 -p 235813 -n tre"
    '''
    call = ["raxml","-T",str(nthreads)," -f d -m GTRCAT -c 25 -p 235813 -n tre -s",aln_file,"> raxml.log"]
    cmd = " ".join(call)
    print("Building a tree via:\n\t" + cmd +
          "\n\tStamatakis, A: RAxML Version 8: A tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies."
          "\n\tIn Bioinformatics, 2014\n")
    os.system(cmd)
    try:
        shutil.copy('RAxML_bestTree.tre', out_file)
        T = Phylo.read(out_file, 'newick')
        if clean_up:
            os.remove('raxml.log')
    except:
        print("TREE BUILDING FAILED, please inspect the raxml.log file\n")
        T=None

    return T

def build_fasttree(aln_file, out_file, clean_up=True):
    '''
    build tree using fasttree with parameters "-nt"
    '''
    call = ["fasttree", "-nt", aln_file, "1>", out_file, "2>", "fasttree.log"]
    cmd = " ".join(call)
    print("Building a tree via:\n\t" + cmd +
          "\n\tPrice et al: FastTree 2 - Approximately Maximum-Likelihood Trees for Large Alignments." +
          "\n\tPLoS ONE 5(3): e9490. https://doi.org/10.1371/journal.pone.0009490\n")
    os.system(cmd)
    try:
        T = Phylo.read(out_file, 'newick')
        if clean_up:
            os.remove('fasttree.log')
    except:
        print("TREE BUILDING FAILED")
        T=None

    return T


def build_iqtree(aln_file, out_file, iqmodel="HKY+F", clean_up=True, nthreads=2):
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

    if iqmodel.lower() != "none":
        call = ["iqtree", "-fast -nt", str(nthreads), "-s", aln_file, "-m", iqmodel,
            ">", "iqtree.log"]
    else:
        call = ["iqtree", "-fast -nt", str(nthreads), "-s", aln_file, ">", "iqtree.log"]

    cmd = " ".join(call)

    print("Building a tree via:\n\t" + cmd +
          "\n\tNguyen et al: IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies."
          "\n\tMol. Biol. Evol., 32:268-274. https://doi.org/10.1093/molbev/msu300\n")
    if iqmodel.lower() == "none":
        print("Conducting a model test... see iqtree.log for the result. You can specify this with --iqmodel in future runs.")
    os.system(cmd)

    # Check result
    try:
        T = Phylo.read(aln_file+".treefile", 'newick')
        shutil.copyfile(aln_file+".treefile", out_file)
        #this allows the user to check intermediate output, as tree.nwk will be
        if clean_up:
            #allow user to see chosen model if modeltest was run
            if iqmodel.lower() == 'none':
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
        with open(stripFile, 'r') as ifile:
            strip_pos = [int(line.strip()) for line in ifile]

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


def run(args):
    # check alignment type, set flags, read in if VCF
    is_vcf = False
    ref = None
    if any([args.alignment.lower().endswith(x) for x in ['.vcf', '.vcf.gz']]):
        if not args.vcf_reference:
            print("ERROR: a reference Fasta is required with VCF-format alignments")
            return -1
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
        variable_fasta = write_out_informative_fasta(compress_seq, args.alignment, stripFile=args.strip_sites)
        fasta = variable_fasta
    else:
        fasta = aln

    if args.iqmodel and not args.method=='iqtree':
        print("Cannot specify model unless using IQTree. Model specification ignored.")

    if args.method=='raxml':
        T = build_raxml(fasta, tree_fname, args.nthreads)
    elif args.method=='iqtree':
        T = build_iqtree(fasta, tree_fname, args.iqmodel, args.nthreads)
    else: #use fasttree - if add more options, put another check here
        T = build_fasttree(fasta, tree_fname)
    end = time.time()
    print("Building original tree took {} seconds".format(str(end-start)))

    if is_vcf and not args.keep_vcf_fasta:
        os.remove(variable_fasta)

    if T:
        import json
        tree_success = Phylo.write(T, tree_fname, 'newick', format_branch_length='%1.8f')
    else:
        return -1
