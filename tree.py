import os, time, sys
sys.path.append('/home/richard/Projects')
from io_util import make_dir, remove_dir, tree_to_json, write_json
from sequences import sequence_set
import numpy as np

class tree(object):
    """tree builds a phylgenetic tree from an alignment and exports it for web visualization"""
    def __init__(self, aln, proteins=None, **kwarks):
        super(tree, self).__init__()
        self.aln = aln
        self.nuc = kwarks['nuc'] if 'nuc' in kwarks else True
        if proteins!=None:
            self.proteins = proteins
        else:
            self.proteins={}
        if 'run_dir' not in kwarks:
            import random
            self.run_dir = '_'.join(['temp', time.strftime('%Y%m%d-%H%M%S',time.gmtime()), str(random.randint(0,1000000))])
        else:
            self.run_dir = kwarks['run_dir']


    def build(self):
        from Bio import Phylo, AlignIO
        from treetime.treetime import io
        from treetime.treetime import utils
        make_dir(self.run_dir)
        os.chdir(self.run_dir)
        for seq in self.aln: seq.name=seq.id
        AlignIO.write(self.aln, 'temp.fasta', 'fasta')

        tree_cmd = ["fasttree"]
        if self.nuc: tree_cmd.append("-nt")
        tree_cmd.append("temp.fasta")
        tree_cmd.append(">")
        tree_cmd.append("initial_tree.nwk")
        os.system(" ".join(tree_cmd))
        self.tt = io.treetime_from_newick('initial_tree.nwk')
        self.tree = self.tt.tree
        io.set_seqs_to_leaves(self.tt, self.aln)
        io.set_node_dates_from_dic(self.tt, {seq.id:utils.numeric_date(seq.attributes['date'])
                                for seq in self.aln if 'date' in seq.attributes})

        os.chdir('..')
        remove_dir(self.run_dir)

    def ancestral(self):
        from treetime.treetime.gtr import GTR
        self.tt.tree.root_at_midpoint()
        self.tt.set_additional_tree_params()
        self.gtr = GTR.standard()
        self.tt.optimize_seq_and_branch_len(self.gtr)


    def timetree(self):
        self.tt.init_date_constraints(self.gtr)
        self.tt.coalescent_model(self.gtr, Tc=0.05)

    def refine(self):
        from treetime.treetime.utils import opt_branch_len
        self.tree.ladderize()
        for node in self.tree.find_clades():
            node.opt_branch_length = opt_branch_len(node)
            if node.up is not None:
                node.muts = ",".join(["".join(map(str, x)) for x in node.mutations])


    def layout(self):
        """Add clade, xvalue, yvalue, mutation and trunk attributes to all nodes in tree"""
        clade = 0
        yvalue = 0
        for node in self.tree.find_clades(order="preorder"):
            node.clade = clade
            clade += 1
            if node.up is not None: #try:
                node.xvalue = node.up.xvalue+node.opt_branch_length
                node.tvalue = node.numdate - self.tree.root.numdate
            else:
                node.xvalue = 0
                node.tvalue = 0
            if node.is_terminal():
                node.yvalue = yvalue
                yvalue += 1
        for node in self.tree.get_nonterminals(order="postorder"):
            node.yvalue = np.mean([x.yvalue for x in node.clades])


    def export(self):
        from Bio import Seq
        from itertools import izip
        timetree_fname = 'tree.json'
        sequence_fname = 'sequences.json'
        tree_json = tree_to_json(self.tree.root)
        write_json(tree_json, timetree_fname, indent=None)
        elems = {}
        elems['root'] = {}
        elems['root']['nuc'] = "".join(self.tree.root.sequence)
        for prot in self.proteins:
            tmp = str(self.proteins[prot].extract(Seq.Seq(elems['root']['nuc'])))
            elems['root'][prot] = str(Seq.translate(tmp.replace('---', 'NNN'))).replace('X','-')


        for node in self.tree.find_clades():
            if hasattr(node, "clade") and hasattr(node, "sequence"):
                elems[node.clade] = {}
                elems[node.clade]['nuc'] = {pos:state for pos, (state, ancstate) in
                                enumerate(izip(node.sequence, self.tree.root.sequence)) if state!=ancstate}
        for prot in self.proteins:
            for node in self.tree.find_clades():
                if hasattr(node, "clade") and hasattr(node, "sequence"):
                    tmp = Seq.translate(str(self.proteins[prot].extract(Seq.Seq("".join(node.sequence)))).replace('-', 'N'))

                    elems[node.clade][prot] = {pos:state for pos, (state, ancstate) in
                                    enumerate(izip(tmp, elems['root'][prot])) if state!=ancstate}

        write_json(elems, sequence_fname, indent=None)


if __name__=="__main__":
    from Bio import SeqIO
    from Bio.SeqFeature import FeatureLocation
    ref_seq = SeqIO.read('NL4-3.gb', 'genbank')
    gene='pol'
    if gene=='gag':
        gag_start = [f.location.start for f in ref_seq.features if f.qualifiers['note'][0]=='gag'][0]
        proteins = {
        'p17': [FeatureLocation(start=f.location.start-gag_start, end=f.location.end-gag_start, strand=1)
                for f in ref_seq.features if f.qualifiers['note'][0]=='p17'][0],
        'p24': [FeatureLocation(start=f.location.start-gag_start, end=f.location.end-gag_start, strand=1)
                for f in ref_seq.features if f.qualifiers['note'][0]=='p24'][0],
        'p6': [FeatureLocation(start=f.location.start-gag_start, end=f.location.end-gag_start, strand=1)
                for f in ref_seq.features if f.qualifiers['note'][0]=='p6'][0],
        'p7': [FeatureLocation(start=f.location.start-gag_start, end=f.location.end-gag_start, strand=1)
                for f in ref_seq.features if f.qualifiers['note'][0]=='p7'][0]}

        myseqs = sequence_set('data/gag.fasta.gz', reference='B|FR|1985|NL4_3_LAI_NY5_pNL43_NL43|244167|NL43|325|U26942')
    elif gene=='pol':
        start = [f.location.start for f in ref_seq.features if f.qualifiers['note'][0]=='pol'][0]
        proteins = {
        'PR': [FeatureLocation(start=f.location.start-start, end=f.location.end-start, strand=1)
              for f in ref_seq.features if f.qualifiers['note'][0]=='PR'][0],
        'RT': [FeatureLocation(start=f.location.start-start, end=f.location.end-start, strand=1)
                for f in ref_seq.features if f.qualifiers['note'][0]=='RT'][0],
        'p15': [FeatureLocation(start=f.location.start-start, end=f.location.end-start, strand=1)
                for f in ref_seq.features if f.qualifiers['note'][0]=='p15'][0],
        'IN': [FeatureLocation(start=f.location.start-start, end=f.location.end-start, strand=1)
                for f in ref_seq.features if f.qualifiers['note'][0]=='IN'][0]}

        myseqs = sequence_set('data/pol.fasta.gz', reference='B|FR|1985|NL4_3_LAI_NY5_pNL43_NL43|244167|NL43|325|U26942')


    myseqs.ungap()
    myseqs.parse({0:"subtype", 1:"country", 2:"date", 4:"name", 5:"id", 6:"patient", 7:"accession"})
    myseqs.parse_date(["%Y-%m-%d", "%Y"])
    myseqs.filter(lambda x:x.attributes['subtype']=='C')
    myseqs.subsample(category = lambda x:x.attributes['date'].year, threshold=10)
    myseqs.codon_align(prune=True)
    myseqs.translate(proteins=proteins)
    myseqs.export_diversity()

    myTree = tree(aln=myseqs.aln, proteins = myseqs.proteins)
    myTree.build()
    myTree.ancestral()
    myTree.timetree()
    myTree.refine()
    myTree.layout()
    myTree.export()

