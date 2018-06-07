import os, shutil, time
from Bio import Phylo
from .utils import read_metadata, get_numerical_dates, write_json
from treetime.vcf_utils import read_vcf, write_vcf
import numpy as np

def timetree(tree=None, aln=None, ref=None, dates=None, keeproot=False, branch_length_mode='auto',
             confidence=False, resolve_polytomies=True, max_iter=2,
             infer_gtr=True, Tc=0.01, reroot='best', use_marginal=False, fixed_pi=None,
             clock_rate=None, n_iqd=None, verbosity=1, **kwarks):
    from treetime import TreeTime

    try: #Tc could be a number or  'opt' or 'skyline'. TreeTime expects a float or int if a number.
        Tc = float(Tc)
    except ValueError:
        True #let it remain a string

    if ref != None: #if VCF, fix pi
        #Otherwise mutation TO gaps is overestimated b/c of seq length
        fixed_pi = [ref.count(base)/len(ref) for base in ['A','C','G','T','-']]
        if fixed_pi[-1] == 0:
            fixed_pi[-1] = 0.05
            fixed_pi = [v-0.01 for v in fixed_pi]

        #set this explicitly if auto, as informative-site only trees can have big branch lengths,
        #making this set incorrectly in TreeTime
        if branch_length_mode == 'auto':
            branch_length_mode = 'joint'

    #send ref, if is None, does no harm
    tt = TreeTime(tree=tree, aln=aln, ref=ref, dates=dates,
                  verbose=verbosity, gtr='JC69')

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



def run(args):
    # check alignment type, set flags, read in if VCF
    is_vcf = False
    ref = None
    tree_meta = {'alignment':args.alignment}
    attributes = ['branch_length']
    # check if tree is provided an can be read
    for fmt in ["newick", "nexus"]:
        try:
            T = Phylo.read(args.tree, fmt)
            tree_meta['input_tree'] = args.tree
            break
        except:
            pass
    if T is None:
        print("ERROR: reading tree from %s failed."%args.tree)
        return -1

    if not args.alignment:
        # fake alignment to appease treetime when only using it for naming nodes...
        if args.ancestral or args.timetree:
            print("ERROR: alignment is required for ancestral reconstruction or timetree inference")
            return -1
        from Bio import SeqRecord, Seq, Align
        seqs = []
        for n in T.get_terminals():
            seqs.append(SeqRecord.SeqRecord(seq=Seq.Seq('ACGT'), id=n.name, name=n.name, description=''))
        aln = Align.MultipleSeqAlignment(seqs)
    elif any([args.alignment.lower().endswith(x) for x in ['.vcf', '.vcf.gz']]):
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


    if args.output:
        tree_fname = args.output
    else:
        tree_fname = '.'.join(args.alignment.split('.')[:-1]) + '_tt.nwk'

    if args.timetree and T:
        if args.metadata is None:
            print("ERROR: meta data with dates is required for time tree reconstruction")
            return -1
        metadata, columns = read_metadata(args.metadata)
        if args.year_limit:
            args.year_limit.sort()
        dates = get_numerical_dates(metadata, fmt=args.date_fmt, min_max_year=args.year_limit)
        for n in T.get_terminals():
            if n.name in metadata and 'date' in metadata[n.name]:
                n.raw_date = metadata[n.name]['date']

        tt = timetree(tree=T, aln=aln, ref=ref, dates=dates, confidence=args.date_confidence,
                      reroot=args.root or 'best',
                      Tc=args.coalescent if args.coalescent is not None else 0.01, #Otherwise can't set to 0
                      use_marginal = args.time_marginal or False,
                      branch_length_mode = args.branch_length_mode or 'auto',
                      clock_rate=args.clock_rate, n_iqd=args.n_iqd)

        tree_meta['clock'] = {'rate':tt.date2dist.clock_rate,
                              'intercept':tt.date2dist.intercept,
                              'rtt_Tmrca':-tt.date2dist.intercept/tt.date2dist.clock_rate}
        attributes.extend(['numdate', 'clock_length', 'mutation_length', 'mutations', 'raw_date', 'date'])
        if not is_vcf:
            attributes.extend(['sequence']) #don't add sequences if VCF - huge!
        if args.date_confidence:
            attributes.append('num_date_confidence')
    elif args.ancestral in ['joint', 'marginal']:
        tt = ancestral_sequence_inference(tree=T, aln=aln, marginal=args.ancestral,
                                          optimize_branch_length=args.branchlengths=='div')
        attributes.extend(['mutation_length', 'mutations', 'sequence'])
    else:
        from treetime import TreeAnc
        # instantiate treetime for the sole reason to name internal nodes
        tt = TreeAnc(tree=T, aln=aln, ref=ref, gtr='JC69', verbose=1)

    if is_vcf:
        #TreeTime overwrites ambig sites on tips during ancestral reconst.
        #Put these back in tip sequences now, to avoid misleading
        tt.recover_var_ambigs()

    tree_meta['nodes'] = prep_tree(T, attributes, is_vcf)

    if T:
        import json
        tree_success = Phylo.write(T, tree_fname, 'newick', format_branch_length='%1.8f')
        if args.node_data:
            node_data_fname = args.node_data
        else:
            node_data_fname = '.'.join(args.alignment.split('.')[:-1]) + '.node_data'

        with open(node_data_fname, 'w') as ofile:
            meta_success = json.dump(tree_meta, ofile)

    #If VCF and ancestral reconst. was done, output VCF including new ancestral seqs
    if is_vcf and (args.ancestral or args.timetree):
        if args.output_vcf:
            vcf_fname = args.output_vcf
        else:
            vcf_fname = '.'.join(args.alignment.split('.')[:-1]) + '.vcf'
        write_vcf(tt.get_tree_dict(keep_var_ambigs=True), vcf_fname)

        return 0 if (tree_success and meta_success) else -1
    else:
        return -1
