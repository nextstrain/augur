import os, shutil, time
from Bio import Phylo
from .utils import read_metadata, get_numerical_dates, write_json

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

    if iqmodel:
        call = ["iqtree", "-fast -nt", str(nthreads), "-s", aln_file, "-m", iqmodel,
            ">", "iqtree.log"]
    else:
        call = ["iqtree", "-fast -nt", str(nthreads), "-s", aln_file, ">", "iqtree.log"]

    cmd = " ".join(call)
    print("Building a tree via:\n\t" + cmd +
          "\n\tNguyen et al: IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies."
          "\n\tMol. Biol. Evol., 32:268-274. https://doi.org/10.1093/molbev/msu300\n")
    os.system(cmd)

    # Check result
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


def timetree(tree=None, aln=None, ref=None, dates=None, keeproot=False,
             confidence=False, resolve_polytomies=True, max_iter=2, dateLimits=None,
             infer_gtr=True, Tc=0.01, reroot='best', use_marginal=False, fixed_pi=None, **kwarks):
    from treetime import TreeTime

    dL_int = None
    if dateLimits:
        dL_int = [int(x) for x in dateLimits]
        dL_int.sort()

    #send ref, if is None, does no harm
    tt = TreeTime(tree=tree, aln=aln, ref=ref, dates=dates,
                  verbose=1, gtr='JC69')

    if confidence and use_marginal:
        # estimate confidence intervals via marginal ML and assign marginal ML times to nodes
        marginal = 'assign'
    else:
        marginal = confidence

    tt.run(infer_gtr=infer_gtr, root=reroot, Tc=Tc, time_marginal=marginal,
           resolve_polytomies=resolve_polytomies, max_iter=max_iter, fixed_pi=fixed_pi, **kwarks)

    if confidence:
        for n in T.find_clades():
            n.numdate_confidence = list(tt.get_max_posterior_region(n, 0.9))

    return tt


def ancestral_sequence_inference(tree=None, aln=None, ref=None, infer_gtr=True,
                                 marginal=False, optimize_branch_length=True):
    from treetime import TreeAnc
    tt = TreeAnc(tree=tree, aln=aln, ref=ref, gtr='JC69', verbose=1)

    if optimize_branch_length:
        tt.optimize_seq_and_branch_len(infer_gtr=infer_gtr, marginal=marginal)
    else: # only infer ancestral sequences, leave branch length untouched
        tt.infer_ancestral_sequences(infer_gtr=infer_gtr, marginal=marginal)

    return tt

def prep_tree(T, attributes):
    data = {}
    for n in T.find_clades():
        data[n.name] = {attr:n.__getattribute__(attr)
                        for attr in attributes if hasattr(n,attr)}
    if 'mutations' in attributes:
        for n in T.find_clades():
            data[n.name]['mutations'] = [[a,int(pos),d] for a,pos,d in data[n.name]['mutations']]
    return data


def run(args):
    # check alignment type, construct reduced alignment if needed
    if any([args.alignment.endswith(x) for x in ['.vcf', '.vcf.gz']]):
        aln = "make alignment from VCF"
    else:
        aln = args.alignment

    T = None
    tree_meta = {'alignment':args.alignment}
    attributes = ['branchlength']
    # check if tree is provided an can be read
    if args.tree:
        for fmt in ["newick", "nexus"]:
            try:
                T = Phylo.read(T, args.tree, fmt)
                tree_meta['input_tree'] = args.tree
                break
            except:
                pass
        if T is None:
            print("WARNING: reading tree from %s failed."
                  "\n\t-- Will attempt to build from alignment."%args.tree)

    start = time.time()
    # without tree, attempt to build tree
    if T is None:
        if args.output:
            tree_fname = args.output
        else:
            tree_fname = '.'.join(args.alignment.split('.')[:-1]) + '.nwk'

        if args.iqmodel and not args.method=='iqtree':
            print("Cannot specify model unless using IQTree. Model specification ignored.")

        tree_meta['topology method'] = args.method
        if args.method=='raxml':
            T = build_raxml(aln, tree_fname, args.nthreads)
        elif args.method=='iqtree':
            T = build_iqtree(aln, tree_fname, args.iqmodel, args.nthreads)
        else: #use fasttree - if add more options, put another check here
            T = build_fasttree(aln, tree_fname)
        end = time.time()
        print("Building original tree took {}".format(str(end-start)))

    if args.timetree and T:
        if args.metadata is None:
            print("ERROR: meta data with dates is required for time tree reconstruction")
            return -1
        metadata, columns = read_metadata(args.metadata)
        dates = get_numerical_dates(metadata, fmt=args.date_fmt)

        tt = timetree(tree=T, aln=aln, dates=dates, confidence=args.date_confidence)
        tree_meta['clock'] = {'rate':tt.date2dist.clock_rate,
                              'intercept':tt.date2dist.intercept,
                              'rtt_Tmrca':-tt.date2dist.intercept/tt.date2dist.clock_rate}
        attributes.extend(['numdate', 'clock_length', 'mutation_length', 'mutations'])
        if args.date_confidence:
            attributes.append('numdate_confidence')
    elif args.ancestral in ['joint', 'marginal']:
        tt = ancestral_sequence_inference(tree=T, aln=aln, marginal=args.ancestral,
                                          optimize_branch_length=args.branchlengths=='div')
        attributes.extend(['mutation_length', 'mutations'])
    else:
        tt = None

    tree_meta['nodes'] = prep_tree(T, attributes)

    if T:
        import json
        tree_success = Phylo.write(T, tree_fname, 'newick')
        if args.node_data:
            node_data_fname = args.node_data
        else:
            node_data_fname = '.'.join(args.alignment.split('.')[:-1]) + '.node_data'

        with open(node_data_fname, 'w') as ofile:
            meta_success = json.dump(tree_meta, ofile)
        return 0 if (tree_success and meta_success) else -1
    else:
        return -1
