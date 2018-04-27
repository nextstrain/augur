import pandas as pd
import numpy as np
import os, shutil, glob, sys
from util import read_sequence_meta_data, parse_date, write_tree_meta_data, write_VCF_style_alignment
from util import collect_tree_meta_data, write_json, generic_argparse, read_in_vcf
from Bio import Phylo


def build_raxml(aln_file, out_file, path, clean_up=True, nthreads=2):
    call = ["raxml","-f d -T",str(nthreads),"-m GTRCAT -c 25 -p 235813 -n tre -s",aln_file,"-w",os.getcwd()+"/"+path+"/results/","> raxml.log"]
    print(" ".join(call))
    os.system(" ".join(call))
    shutil.copy(path+'/results/RAxML_bestTree.tre', out_file)
    try:
        T = Phylo.read(out_file, 'newick')
        if clean_up:
            os.remove('raxml.log')
    except:
        print("TREE BUILDING FAILED")
        T=None

    return T

def build_fasttree(aln_file, out_file, clean_up=True):
    call = ["fasttree", "-nt", aln_file, "1>", out_file, "2>", "fasttree.log"]
    print(" ".join(call))
    os.system(" ".join(call))
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
    if iqmodel:
        call = ["iqtree", "-nt", str(nthreads), "-s", aln_file, "-m", iqmodel[0],
            ">", "iqtree.log"]
    else:
        call = ["iqtree", "-nt", str(nthreads), "-s", aln_file, ">", "iqtree.log"]
    print(" ".join(call))
    os.system(" ".join(call))
    try:
        T = Phylo.read(aln_file+".treefile", 'newick')
        shutil.copyfile(aln_file+".treefile", out_file)
        #this allows the user to check intermediate output, as tree.nwk will be
        #written over with TreeTime tree
        shutil.copyfile(aln_file+".treefile", out_file.replace(".nwk",".iqtree.nwk"))
        if clean_up:
            #allow user to see chosen model
            shutil.copyfile('iqtree.log', out_file.replace("tree.nwk","iqtree.log"))
            os.remove('iqtree.log')
            for filename in glob.glob(aln_file+".*"):
                os.remove(filename)
    except:
        print("TREE BUILDING FAILED")
        T=None
    return T


def timetree(tree=None, aln=None, ref=None, seq_meta=None, keeproot=False,
             confidence=False, resolve_polytomies=True, max_iter=2, dateLimits=None,
             infer_gtr=True, Tc=0.01, reroot='best', use_marginal=False, **kwarks):
    from treetime import TreeTime

    dL_int = None
    if dateLimits:
        dL_int = [int(x) for x in dateLimits]
        dL_int.sort()

    dates = {}
    for name, data in seq_meta.items():
        num_date = parse_date(data["date"], date_fmt, dL_int)
        if num_date is not None:
            dates[name] = num_date

    #send ref, if is None, does no harm
    tt = TreeTime(tree=tree, aln=aln, ref=ref, dates=dates, gtr='JC69')

    if confidence and use_marginal:
        # estimate confidence intervals via marginal ML and assign marginal ML times to nodes
        marginal = 'assign'
    else:
        marginal = confidence

    #Length of VCF files means GTR model with gaps causes overestimation of mutation TO gaps
    #so gaps appear in internal nodes when no gaps at tips! To prevent....
    pi = None
    if ref != None: #if VCF, fix pi
        pi = np.array([0.1618, 0.3188, 0.3176, 0.1618, 0.04]) #from real runs (Walker 2018)


    tt.run(infer_gtr=infer_gtr, root=reroot, Tc=Tc, time_marginal=marginal,
           resolve_polytomies=resolve_polytomies, max_iter=max_iter, fixed_pi=pi, **kwarks)





    for n in T.find_clades():
        n.num_date = n.numdate # treetime convention is different from augur...
        # get 90% max posterior region)
        if confidence:
            n.num_date_confidence = list(tt.get_max_posterior_region(n, 0.9))
    return tt


def ancestral_sequence_inference(tree=None, aln=None, ref=None, infer_gtr=True,
                                 optimize_branch_length=True):
    from treetime import TreeAnc
    tt = TreeAnc(tree=tree, aln=aln, ref=ref, gtr='JC69')

    if optimize_branch_length:
        tt.optimize_seq_and_branch_len(infer_gtr=infer_gtr)
    else: # only infer ancestral sequences, leave branch length untouched
        tt.infer_ancestral_sequences(infer_gtr=infer_gtr)

    return tt


def export_sequence_fasta(T, path):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
    from Bio import AlignIO

    fname = tree_sequence_alignment(path, 'nuc')
    seqs = [SeqRecord(Seq(''.join(T.root.sequence)), name='root', id='root')]
    for node in T.find_clades():
        seqs.append(SeqRecord(Seq(''.join(node.sequence)), description='', name=node.name, id=node.name))
    AlignIO.write(MultipleSeqAlignment(seqs), fname, 'fasta')



def main()
    parser = generic_argparse("Build the tree from the prepared sequence data")
    parser.add_argument('-s', help='alignment')
    parser.add_argument('--tree_out', help='output file')
    parser.add_argument('--tree_meta', help='output file')
    parser.add_argument('--nthreads', type=int, default=2,
                        help='number of threads')
    parser.add_argument('--ancestral', action='store_true', default=False,
                        help='calculate and store ancestral sequences')
    parser.add_argument('--timetree', action='store_true', default=False,
                       help='infer time stamped phylogeny')
    parser.add_argument('--confidence', action='store_true', default=False,
                       help='estimate confidence intervals for node timing')
    parser.add_argument('--Tc', type=float, default=0.0,
                       help='coalescence time scale measured in substitution rate units')
    parser.add_argument('--keeproot', action='store_true', default=False,
                        help="don't reroot the tree")

    #EBH 4 Dec 2017
    parser.add_argument('--iqtree', action='store_true', default=False,
                        help="use iqtree for initial tree building")
    parser.add_argument('--raxml', action='store_true', default=False,
                        help="use raxml for initial tree building")
    parser.add_argument('--fasttree', action='store_true', default=True,
                        help="use fasttree for initial tree building (default)")
    parser.add_argument('--vcf', action='store_true', default=False,
                        help="sequence is in VCF format")

    #EBH 5 Jan 2018
    parser.add_argument('--iqmodel', nargs=1, help='model to use with iqtree')
    parser.add_argument('--drm', type=str,
                        help="file of DRMs to exclude from inital tree-building")

    #EBH 14 Feb 2018
    parser.add_argument('--roottype', nargs="+",#type=str, default="residual",
                        help="type of rerooting. options are 'rsq', 'residual' (default), and 'oldest'")

    #EBH 16 Mar 2018
    parser.add_argument('--varAmbigs', action='store_true', default=False,
                        help="preserve ambiguous bases at variable sites in recorded mutations")
    parser.add_argument('--dateLimit', nargs='+',
                        help="specify min and max year for samples without dates. Order doesn't matter. If only one value, taken as min, and max set to current year.")

    args = parser.parse_args()
    tree_fname = args.tree_out
    tree_prop_fname = args.tree_meta


    date_fmt = '%Y-%m-%d'

    import time

    if args.iqmodel and not args.iqtree:
        print("Cannot specify model unless using IQTree. Model specification ignored.")

    start = time.time()
    if args.raxml:
        T = build_raxml(args.s, tree_fname, path, args.nthreads)
    elif args.iqtree:
        T = build_iqtree(args.s, tree_fname, args.iqmodel, args.nthreads)
    else: #use fasttree - if add more options, put another check here
        T = build_fasttree(args.s, tree_fname)
    end = time.time()
    print "Building original tree took {}".format(str(end-start))

    meta = read_sequence_meta_data(path)
    fields = ['branchlength', 'clade']

    #Anything but a list of sequences to root by, shouldn't go as a "list".
    if len(args.roottype) == 1:
        args.roottype = args.roottype[0]

    start = time.time()
    if args.timetree:
        if args.vcf:
            tt = timetree(tree=T, aln=sequences, ref=ref, confidence=args.confidence, dateLimits=args.dateLimit,
                          seq_meta=meta, reroot=None if args.keeproot else args.roottype, Tc=args.Tc)#, use_marginal=True)
        else:
            tt = timetree(tree=T, aln=ref_alignment(path), confidence=args.confidence, dateLimits=args.dateLimit,
                          seq_meta=meta, reroot=None if args.keeproot else args.roottype, Tc=args.Tc)

        T = tt.tree
        fields.extend(['mutations', 'mutation_length', 'num_date', 'clock_length'])
        if args.confidence:
            fields.append('num_date_confidence')
    elif args.ancestral:
        if args.vcf:
            tt = ancestral_sequence_inference(tree=T, aln=sequences, ref=ref)
        else:
            tt = ancestral_sequence_inference(tree=T, aln=ref_alignment(path))
        T = tt.tree
        fields.extend(['mutations', 'mutation_length'])

    end = time.time()
    print "TreeTime took {}".format(str(end-start))

    clade_index = 0
    for n in T.find_clades(order='preorder'):
        n.clade = clade_index
        clade_index+=1

    if args.vcf:
        Phylo.write(T, tree_fname, 'newick', format_branch_length='%1.8f')
        #This gives more digits to branch length which for small numbers
        #means it's more than just 0 and 0.00001 in the Newick!
    else:
        Phylo.write(T, tree_fname, 'newick')

    if args.varAmbigs: #if requested, put ambigs back on tips
        tt.recover_var_ambigs()
    meta_dic = collect_tree_meta_data(T, fields, args.vcf)
    write_tree_meta_data(path, meta_dic)

    with open(sequence_gtr_model(path),'w') as ofile:
        ofile.write(str(tt.gtr))

    start = time.time()
    #do NOT print out all full sequences if VCF - will be huge!
    if args.timetree or args.ancestral:
        if args.vcf:
            export_sequence_VCF(tt, path, args.varAmbigs)
        else:
            export_sequence_fasta(T, path)
    end = time.time()
    print "Writing out VCF/Fasta took {}".format(str(end-start))


if __name__ == '__main__':
    main(sys.argv)