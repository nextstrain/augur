from __future__ import division, print_function
import os, time
from io_util import make_dir, remove_dir, tree_to_json, write_json
from sequences import sequence_set
import numpy as np

def resolve_polytomies(tree):
    for node in tree.get_nonterminals('preorder'):
        node.confidence = None
        if len(node.clades)>2:
            n = len(node.clades)
            children = list(node.clades)
            node.clades = []
            node.split(branch_length=1e-5)
            if n>3:
                node.clades[0].clades = children[:len(children)//2]
                node.clades[1].clades = children[len(children)//2:]
                for c in node.clades:
                    c.name=''
                    c.confidence = None
            else:
                node.clades[0] = children[0]
                node.clades[1].clades = children[1:]
                node.clades[1].confidence = None
                node.clades[1].name = None


class tree(object):
    """tree builds a phylgenetic tree from an alignment and exports it for web visualization"""
    def __init__(self, aln, proteins=None, **kwarks):
        super(tree, self).__init__()
        self.aln = aln
        self.nthreads = 2
        self.sequence_lookup = {seq.id:seq for seq in aln}
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


    def build(self, root='midpoint', raxml=True, raxml_time_limit=0.5):
        from Bio import Phylo, AlignIO
        import subprocess, glob, shutil
        make_dir(self.run_dir)
        os.chdir(self.run_dir)
        for seq in self.aln: seq.name=seq.id
        AlignIO.write(self.aln, 'temp.fasta', 'fasta')

        tree_cmd = ["fasttree"]
        if self.nuc: tree_cmd.append("-nt")
        tree_cmd.append("temp.fasta")
        tree_cmd.append(">")
        tree_cmd.append("initial_tree.newick")
        os.system(" ".join(tree_cmd))

        out_fname = "tree_infer.newick"
        if raxml:
            if raxml_time_limit>0:
                tmp_tree = Phylo.read('initial_tree.newick','newick')
                resolve_iter = 0
                resolve_polytomies(tmp_tree)
                while (not tmp_tree.is_bifurcating()) and (resolve_iter<10):
                    resolve_iter+=1
                    resolve_polytomies(tmp_tree)
                Phylo.write(tmp_tree,'initial_tree.newick', 'newick')
                AlignIO.write(self.aln,"temp.phyx", "phylip-relaxed")
                print( "RAxML tree optimization with time limit", raxml_time_limit,  "hours")
                # using exec to be able to kill process
                end_time = time.time() + int(raxml_time_limit*3600)
                process = subprocess.Popen("exec raxml -f d -T " + str(self.nthreads) + " -j -s temp.phyx -n topology -c 25 -m GTRCAT -p 344312987 -t initial_tree.newick", shell=True)
                while (time.time() < end_time):
                    if os.path.isfile('RAxML_result.topology'):
                        break
                    time.sleep(10)
                process.terminate()

                checkpoint_files = glob.glob("RAxML_checkpoint*")
                if os.path.isfile('RAxML_result.topology'):
                    checkpoint_files.append('RAxML_result.topology')
                if len(checkpoint_files) > 0:
                    last_tree_file = checkpoint_files[-1]
                    shutil.copy(last_tree_file, 'raxml_tree.newick')
                else:
                    shutil.copy("initial_tree.newick", 'raxml_tree.newick')
            else:
                shutil.copy("initial_tree.newick", 'raxml_tree.newick')

            try:
                print("RAxML branch length optimization")
                os.system("raxml -f e -T " + str(self.nthreads) + " -s temp.phyx -n branches -c 25 -m GTRGAMMA -p 344312987 -t raxml_tree.newick")
                shutil.copy('RAxML_result.branches', out_fname)
            except:
                print("RAxML branch length optimization failed")
                shutil.copy('raxml_tree.newick', out_fname)
        else:
            shutil.copy('initial_tree.newick', out_fname)
        self.tt_from_file(out_fname, root)
        os.chdir('..')
        remove_dir(self.run_dir)
        self.is_timetree=False


    def tt_from_file(self, infile, root='none'):
        from treetime.gtr import GTR
        from treetime import io, utils
        gtr = GTR.standard()
        self.tt = io.treetime_from_newick(gtr, infile)
        io.set_seqs_to_leaves(self.tt, self.aln)
        io.set_node_dates_from_dic(self.tt, {seq.id:utils.numeric_date(seq.attributes['date'])
                                for seq in self.aln if 'date' in seq.attributes})
        self.tree = self.tt.tree
        if root=='midpoint':
            self.tt.tree.root_at_midpoint()
            self.tt.set_additional_tree_params()
        elif root=='oldest':
            tmp = self.tt.reroot_to_oldest()

        for node in self.tree.get_terminals():
            if node.name in self.sequence_lookup:
                seq = self.sequence_lookup[node.name]
                for attr in seq.attributes:
                    if attr == 'date':
                        node.date = seq.attributes['date'].strftime('%Y-%m-%d')
                    else:
                        node.__setattr__(attr, seq.attributes[attr])


    def ancestral(self):
        self.tt.optimize_seq_and_branch_len(infer_gtr=True)

    def timetree(self, Tc=0.05, infer_gtr=False, **kwarks):
        print('rerooting...')
        self.tt.reroot_to_best_root(infer_gtr=True, **kwarks)
        print('estimating time tree with coalescent model...')
        self.tt.coalescent_model(Tc=Tc,**kwarks)

        self.is_timetree=True

    def geo_inference(self, attr):
        from treetime.gtr import GTR
        places = set()
        for node in self.tree.find_clades():
            if hasattr(node, attr):
                places.add(node.__getattribute__(attr))
            if hasattr(node, 'sequence'):
                node.nuc_sequence = node.sequence
        places = sorted(places)
        alphabet = {chr(65+i):place for i,place in enumerate(places)}
        alphabet_rev = {v:k for k,v in alphabet.iteritems()}
        self.tt.sequence_gtr = self.tt.gtr
        nc = len(places)
        myGeoGTR = GTR.custom(pi = np.ones(nc, dtype=float)/nc, W=np.ones((nc,nc)),
                              alphabet = np.array(sorted(alphabet.keys())))
        myGeoGTR.profile_map['-'] = np.ones(nc)
        for node in self.tree.get_terminals():
            if hasattr(node, attr):
                node.sequence=np.array([alphabet_rev[node.__getattribute__(attr)]])
            else:
                node.sequence=np.array(['-'])
        for node in self.tree.get_nonterminals():
            node.__delattr__('sequence')
        self.tt._gtr = myGeoGTR
        self.tt.reconstruct_anc(method='ml')
        self.tt.infer_gtr(print_raw=True)
        self.tt.reconstruct_anc(method='ml')

        self.tt.geogtr = self.tt.gtr
        self.tt.geogtr.alphabet_to_location = alphabet
        self.tt._gtr = self.tt.sequence_gtr
        for node in self.tree.find_clades():
            node.__setattr__(attr, alphabet[node.sequence[0]])
            node.sequence = node.nuc_sequence


    def add_translations(self):
        from Bio import Seq
        for node in self.tree.find_clades():
            if not hasattr(node, "translations"):
                node.translations={}
            for prot in self.proteins:
                node.translations[prot] = Seq.translate(str(self.proteins[prot].extract(Seq.Seq("".join(node.sequence)))).replace('-', 'N'))

    def refine(self):
        from treetime.utils import opt_branch_len
        self.tree.ladderize()
        for node in self.tree.find_clades():
            node.opt_branch_length = opt_branch_len(node)
            if node.up is not None:
                node.mut_str = ",".join(["".join(map(str, x)) for x in node.mutations])
                node.aa_mut_str = {}
                node.aa_mutations = {}
                if hasattr(node, 'translations'):
                    for prot in node.translations:
                        node.aa_mutations[prot] = [(a,pos+1,d) for pos, (a,d) in
                                            enumerate(zip(node.up.translations[prot],
                                                          node.translations[prot])) if a!=d]
                        node.aa_mut_str[prot] = ",".join(["".join(map(str,x)) for x in node.aa_mutations[prot]])

    def layout(self):
        """Add clade, xvalue, yvalue, mutation and trunk attributes to all nodes in tree"""
        clade = 0
        yvalue = self.tree.count_terminals()
        for node in self.tree.find_clades(order="preorder"):
            node.clade = clade
            clade += 1
            if node.up is not None: #try:
                if self.is_timetree:
                    node.xvalue = node.up.xvalue+node.opt_branch_length
                    node.tvalue = node.numdate - self.tree.root.numdate
                else:
                    node.xvalue = node.up.xvalue+node.branch_length
                    node.tvalue = 0
            else:
                node.xvalue = 0
                node.tvalue = 0
            if node.is_terminal():
                node.yvalue = yvalue
                yvalue -= 1
        for node in self.tree.get_nonterminals(order="postorder"):
            node.yvalue = np.mean([x.yvalue for x in node.clades])


    def export(self, path = '', extra_attr = ['aa_mut_str']):
        from Bio import Seq
        from itertools import izip
        timetree_fname = path+'tree.json'
        sequence_fname = path+'sequences.json'
        tree_json = tree_to_json(self.tree.root, extra_attr=extra_attr)
        write_json(tree_json, timetree_fname, indent=None)
        elems = {}
        elems['root'] = {}
        elems['root']['nuc'] = "".join(self.tree.root.sequence)
        for prot in self.proteins:
            tmp = str(self.proteins[prot].extract(Seq.Seq(elems['root']['nuc'])))
            #elems['root'][prot] = str(Seq.translate(tmp.replace('---', 'NNN'))).replace('X','-')
            elems['root'][prot] = str(Seq.translate(tmp.replace('-', 'N'))).replace('X','-')


        for node in self.tree.find_clades():
            if hasattr(node, "clade") and hasattr(node, "sequence"):
                elems[node.clade] = {}
                elems[node.clade]['nuc'] = {pos:state for pos, (state, ancstate) in
                                enumerate(izip(node.sequence, self.tree.root.sequence)) if state!=ancstate}
        for node in self.tree.find_clades():
            if hasattr(node, "clade") and hasattr(node, "translations"):
                for prot in self.proteins:
                    elems[node.clade][prot] = {pos:state for pos, (state, ancstate) in
                                    enumerate(izip(node.translations[prot], elems['root'][prot])) if state!=ancstate}

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

