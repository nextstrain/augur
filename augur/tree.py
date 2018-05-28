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


def timetree(tree=None, aln=None, ref=None, dates=None, keeproot=False, branch_length_mode='auto',
             confidence=False, resolve_polytomies=True, max_iter=2, 
             infer_gtr=True, Tc=0.01, reroot='best', use_marginal=False, fixed_pi=None,
             clock_rate=None, n_iqd=None, **kwarks):
    from treetime import TreeTime

    if ref != None: #if VCF, fix pi
        #Otherwise mutation TO gaps is overestimated b/c of seq length
        fixed_pi = [ref.count(base)/len(ref) for base in ['A','C','G','T','-']]
        if fixed_pi[-1] == 0:
            fixed_pi[-1] = 0.05
            fixed_pi = [v-0.01 for v in fixed_pi]
        
        #set this explicitly, as informative-site only trees can have big branch lengths,
        #making this set incorrectly in TreeTime
        branch_length_mode = 'joint'

    #send ref, if is None, does no harm
    tt = TreeTime(tree=tree, aln=aln, ref=ref, dates=dates,
                  verbose=1, gtr='JC69')

    if confidence and use_marginal:
        # estimate confidence intervals via marginal ML and assign marginal ML times to nodes
        marginal = 'assign'
    else:
        marginal = confidence

    tt.run(infer_gtr=infer_gtr, root=reroot, Tc=Tc, time_marginal=marginal,
           branch_length_mode=branch_length_mode, resolve_polytomies=resolve_polytomies,
           max_iter=max_iter, fixed_pi=fixed_pi, fixed_clock_rate=clock_rate,
           n_iqd=n_iqd, **kwarks)

    if confidence:
        for n in tt.tree.find_clades():
            n.num_date_confidence = list(tt.get_max_posterior_region(n, 0.9))

    print("\nInferred a time resolved phylogeny using TreeTime:"
          "\n\tSagulenko et al. TreeTime: Maximum-likelihood phylodynamic analysis"
          "\n\tVirus Evolution, vol 4, https://academic.oup.com/ve/article/4/1/vex042/4794731\n")
    return tt


def ancestral_sequence_inference(tree=None, aln=None, ref=None, infer_gtr=True,
                                 marginal=False, optimize_branch_length=True):
    from treetime import TreeAnc
    tt = TreeAnc(tree=tree, aln=aln, ref=ref, gtr='JC69', verbose=1)

    if optimize_branch_length:
        tt.optimize_seq_and_branch_len(infer_gtr=infer_gtr, marginal=marginal)
    else: # only infer ancestral sequences, leave branch length untouched
        tt.infer_ancestral_sequences(infer_gtr=infer_gtr, marginal=marginal)

    print("\nInferred ancestral sequence states using TreeTime:"
          "\n\tSagulenko et al. TreeTime: Maximum-likelihood phylodynamic analysis"
          "\n\tVirus Evolution, vol 4, https://academic.oup.com/ve/article/4/1/vex042/4794731\n")

    return tt

def prep_tree(T, attributes, is_vcf=False):
    data = {}
    inc = 1 if is_vcf else 0 #convert python numbering to start-at-1
    for n in T.find_clades():
        data[n.name] = {attr:n.__getattribute__(attr)
                        for attr in attributes if hasattr(n,attr)}
    if 'mutations' in attributes:
        for n in T.find_clades():
            data[n.name]['mutations'] = [[a,int(pos)+inc,d] for a,pos,d in data[n.name]['mutations']]
    if not is_vcf and 'sequence' in attributes: #don't attach sequence if VCF!
        for n in T.find_clades():
            if hasattr(n, 'sequence'):
                data[n.name]['sequence'] = ''.join(n.sequence)
            else:
                data[n.name]['sequence']=''

    return data

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

    T = None
    tree_meta = {'alignment':args.alignment}
    attributes = ['branchlength']
    # check if tree is provided an can be read
    if args.tree:
        for fmt in ["newick", "nexus"]:
            try:
                T = Phylo.read(args.tree, fmt)
                tree_meta['input_tree'] = args.tree
                break
            except:
                pass
        if T is None:
            print("WARNING: reading tree from %s failed."
                  "\n\t-- Will attempt to build from alignment."%args.tree)

    start = time.time()

    if args.output:
        tree_fname = args.output
    else:
        tree_fname = '.'.join(args.alignment.split('.')[:-1]) + '.nwk'

    # without tree, attempt to build tree
    if T is None:
        # construct reduced alignment if needed
        if is_vcf:
            variable_fasta = write_out_informative_fasta(compress_seq, args.alignment, stripFile=args.strip_sites)
            fasta = variable_fasta
        else:
            fasta = aln

        if args.iqmodel and not args.method=='iqtree':
            print("Cannot specify model unless using IQTree. Model specification ignored.")

        tree_meta['topology method'] = args.method
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

    if args.timetree and T:
        if args.metadata is None:
            print("ERROR: meta data with dates is required for time tree reconstruction")
            return -1
        metadata, columns = read_metadata(args.metadata)
        if args.year_limit: 
            args.year_limit.sort()
        dates = get_numerical_dates(metadata, fmt=args.date_fmt, min_max_year=args.year_limit)

        tt = timetree(tree=T, aln=aln, ref=ref, dates=dates, confidence=args.date_confidence,
                      reroot=args.root if args.root else 'best', 
                      clock_rate=args.clock_rate, n_iqd=args.n_iqd)

        tree_meta['clock'] = {'rate':tt.date2dist.clock_rate,
                              'intercept':tt.date2dist.intercept,
                              'rtt_Tmrca':-tt.date2dist.intercept/tt.date2dist.clock_rate}
        attributes.extend(['numdate', 'clock_length', 'mutation_length', 'mutations'])
        if not is_vcf:
            attributes.extend(['sequence']) #don't add sequences if VCF - huge!
        if args.date_confidence:
            attributes.append('num_date_confidence')
    elif args.ancestral in ['joint', 'marginal']:
        tt = ancestral_sequence_inference(tree=T, aln=aln, marginal=args.ancestral,
                                          optimize_branch_length=args.branchlengths=='div')
        attributes.extend(['mutation_length', 'mutations', 'sequence'])
    else:
        tt = None

    if is_vcf:
        #TreeTime overwrites ambig sites on tips during ancestral reconst.
        #Put these back in tip sequences now, to avoid misleading
        tt.recover_var_ambigs()

    tree_meta['nodes'] = prep_tree(T, attributes, is_vcf)

    if T:
        import json
        tree_success = Phylo.write(T, tree_fname, 'newick', format_branch_length='%1.8f')
        if args.timetree or args.ancestral in ['joint', 'marginal']:
            if args.node_data:
                node_data_fname = args.node_data
            else:
                node_data_fname = '.'.join(args.alignment.split('.')[:-1]) + '.node_data'

            with open(node_data_fname, 'w') as ofile:
                meta_success = json.dump(tree_meta, ofile)
        else:
            meta_success=True

    #If VCF and ancestral reconst. was done, output VCF including new ancestral seqs
    if is_vcf and (args.ancestral or args.treetime):
        if args.output_vcf:
            vcf_fname = args.output_vcf
        else:
            vcf_fname = '.'.join(args.alignment.split('.')[:-1]) + '.vcf'
        write_vcf(tt.get_tree_dict(keep_var_ambigs=True), vcf_fname)

        return 0 if (tree_success and meta_success) else -1
    else:
        return -1
