import os, shutil, time
from Bio import Phylo

def build_raxml(aln_file, out_file, clean_up=True, nthreads=2):
    call = ["raxml","-f d -T",str(nthreads),"-m GTRCAT -c 25 -p 235813 -n tre -s",aln_file,"> raxml.log"]
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
    call = ["fasttree", "-nt", aln_file, "1>", out_file, "2>", "fasttree.log"]
    cmd = " ".join(call)
    print("Building a tree via:\n\t" + cmd +
          "\n\tPrice et al: FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments."
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


def build_iqtree(aln_file, out_file, iqmodel, clean_up=True, nthreads=2):
    #return Phylo.read(out_file.replace(".nwk",".iqtree.nwk"), 'newick') #uncomment for debug skip straight to TreeTime

    with open(aln_file) as ifile:
        tmp_seqs = ifile.readlines()

    aln_file = "temp_iqtree.fasta"
    with open(aln_file, 'w') as ofile:
        for line in tmp_seqs:
            ofile.write(line.replace('/', '_X_X_').replace('|','_Y_Y_'))

    if iqmodel:
        call = ["iqtree", "-nt", str(nthreads), "-s", aln_file, "-m", iqmodel[0],
            ">", "iqtree.log"]
    else:
        call = ["iqtree", "-fast -nt", str(nthreads), "-s", aln_file, ">", "iqtree.log"]
    cmd = " ".join(call)
    print("Building a tree via:\n\t" + cmd +
          "\n\tNguyen et al: IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies."
          "\n\tMol. Biol. Evol., 32:268-274. https://doi.org/10.1093/molbev/msu300\n")
    os.system(cmd)
    try:
        T = Phylo.read(aln_file+".treefile", 'newick')
        shutil.copyfile(aln_file+".treefile", out_file)
        #this allows the user to check intermediate output, as tree.nwk will be
        if clean_up:
            #allow user to see chosen model
            shutil.copyfile('iqtree.log', out_file.replace("tree.nwk","iqtree.log"))
            os.remove('iqtree.log')
            os.remove(aln_file)
            for ext in [".bionj",".ckp.gz",".iqtree",".log",".mldist",".model.gz",".treefile"]:
                if os.path.isfile(aln_file + ext):
                    os.remove(aln_file + ext)
            for n in T.get_terminals():
                n.name = n.name.replace('_X_X_','/').replace('_Y_Y_','|')
    except:
        print("TREE BUILDING FAILED")
        T=None
    return T

def run(args):
    if args.output:
        tree_fname = args.output
    else:
        tree_fname = '.'.join(args.alignment.split('.')) + '.nwk'
    aln = args.alignment

    if args.iqmodel and not args.method=='iqtree':
        print("Cannot specify model unless using IQTree. Model specification ignored.")

    start = time.time()
    if args.method=='raxml':
        T = build_raxml(aln, tree_fname, args.nthreads)
    elif args.method=='iqtree':
        T = build_iqtree(aln, tree_fname, args.iqmodel, args.nthreads)
    else: #use fasttree - if add more options, put another check here
        T = build_fasttree(aln, tree_fname)
    end = time.time()
    print("Building original tree took {}".format(str(end-start)))
    if T:
        return Phylo.write(T, tree_fname, 'newick')
    else:
        return 0