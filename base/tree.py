from __future__ import division, print_function
import os, time
from io_util import make_dir, remove_dir, tree_to_json, write_json, myopen
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
        self.dump_attr = []
        if proteins!=None:
            self.proteins = proteins
        else:
            self.proteins={}
        if 'run_dir' not in kwarks:
            import random
            self.run_dir = '_'.join(['temp', time.strftime('%Y%m%d-%H%M%S',time.gmtime()), str(random.randint(0,1000000))])
        else:
            self.run_dir = kwarks['run_dir']


    def dump(self, treefile, nodefile):
        from Bio import Phylo
        Phylo.write(self.tree, treefile, 'newick')
        node_props = {}
        for node in self.tree.find_clades():
            node_props[node.name] = {attr:node.__getattribute__(attr) for attr in self.dump_attr if hasattr(node, attr)}

        with myopen(nodefile, 'w') as nfile:
            from cPickle import dump
            dump(node_props, nfile)


    def build(self, root='midpoint', raxml=True, raxml_time_limit=0.5, raxml_bin='raxml'):
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
                process = subprocess.Popen("exec " + raxml_bin + " -f d -T " + str(self.nthreads) + " -j -s temp.phyx -n topology -c 25 -m GTRCAT -p 344312987 -t initial_tree.newick", shell=True)
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
                os.system(raxml_bin + " -f e -T " + str(self.nthreads) + " -s temp.phyx -n branches -c 25 -m GTRGAMMA -p 344312987 -t raxml_tree.newick")
                shutil.copy('RAxML_result.branches', out_fname)
            except:
                print("RAxML branch length optimization failed")
                shutil.copy('raxml_tree.newick', out_fname)
        else:
            shutil.copy('initial_tree.newick', out_fname)
        self.tt_from_file(out_fname, root)
        os.chdir('..')
        remove_dir(self.run_dir)


    def tt_from_file(self, infile, root='best', nodefile=None):
        from treetime import TreeTime
        from treetime import utils
        self.is_timetree=False
        print('Reading tree from file',infile)
        dates  =   {seq.id:seq.attributes['num_date']
                    for seq in self.aln if 'date' in seq.attributes}
        self.tt = TreeTime(dates=dates, tree=infile, gtr='Jukes-Cantor', aln = self.aln, verbose=4)
        if root:
            self.tt.reroot(root=root)
        self.tree = self.tt.tree

        for node in self.tree.find_clades():
            if node.is_terminal() and node.name in self.sequence_lookup:
                seq = self.sequence_lookup[node.name]
                node.attr = seq.attributes
                try:
                    node.attr['date'] = node.attr['date'].strftime('%Y-%m-%d')
                except:
                    pass
            else:
                node.attr = {}

        if nodefile is not None:
            print('reading node properties from file:',nodefile)
            with myopen(nodefile, 'r') as infile:
                from cPickle import load
                node_props = load(infile)
            for n in self.tree.find_clades():
                if n.name in node_props:
                    for attr in node_props[n.name]:
                        n.__setattr__(attr, node_props[n.name][attr])
                else:
                    print("No node properties found for ", n.name)


    def ancestral(self):
        self.tt.optimize_seq_and_branch_len(infer_gtr=True)
        self.dump_attr.append('sequence')


    def timetree(self, Tc=0.01, infer_gtr=True, reroot='best', resolve_polytomies=True, max_iter=2, **kwarks):
        self.tt.run(infer_gtr=infer_gtr, root=reroot, Tc=Tc,
                    resolve_polytomies=resolve_polytomies, max_iter=max_iter)
        print('estimating time tree...')
        self.dump_attr.extend(['numdate','date','sequence'])
        for node in self.tree.find_clades():
            if hasattr(node,'attr'):
                node.attr['num_date'] = node.numdate
            else:
                node.attr = {'num_date':node.numdate}
        self.is_timetree=True


    def geo_inference(self, attr):
        '''
        infer a "mugration" model by pretending each region corresponds to a sequence
        state and repurposing the GTR inference and ancestral reconstruction
        '''
        from treetime.gtr import GTR
        # Determine alphabet and store reconstructed ancestral sequences
        places = set()
        nuc_seqs = {}
        nuc_muts = {}
        for node in self.tree.find_clades():
            if hasattr(node, 'attr'):
                if attr in node.attr:
                    places.add(node.attr[attr])
            if hasattr(node, 'sequence'):
                nuc_seqs[node] = node.sequence
            if hasattr(node, 'mutations'):
                nuc_muts[node] = node.mutations
                node.__delattr__('mutations')

        # construct GTR (flat for now). The missing DATA symbol is a '-' (ord('-')==45)
        places = sorted(places)
        nc = len(places)
        if nc<2 or nc>180:
            print("geo_inference: can't have less than 2 or more than 180 places!")
            return
        alphabet = {chr(65+i):place for i,place in enumerate(places)}
        alphabet_rev = {v:k for k,v in alphabet.iteritems()}
        sequence_gtr = self.tt.gtr
        myGeoGTR = GTR.custom(pi = np.ones(nc, dtype=float)/nc, W=np.ones((nc,nc)),
                              alphabet = np.array(sorted(alphabet.keys())))
        myGeoGTR.profile_map['-'] = np.ones(nc)

        # set geo info to nodes as one letter sequence.
        for node in self.tree.get_terminals():
            if hasattr(node, 'attr'):
                if attr in node.attr:
                    node.sequence=np.array([alphabet_rev[node.attr[attr]]])
            else:
                node.sequence=np.array(['-'])
        for node in self.tree.get_nonterminals():
            node.__delattr__('sequence')
        # set custom GTR model, run inference
        self.tt._gtr = myGeoGTR
        tmp_use_mutation_length = self.tt.use_mutation_length
        self.tt.use_mutation_length=False
        self.tt.infer_ancestral_sequences(method='ml', infer_gtr=True,
            store_compressed=False, pc=5.0, marginal=True)

        # restore the nucleotide sequence and mutations to maintain expected behavior
        self.tt.geogtr = self.tt.gtr
        self.tt.geogtr.alphabet_to_location = alphabet
        self.tt._gtr = sequence_gtr
        self.dump_attr.append(attr)
        for node in self.tree.find_clades():
            node.attr[attr] = alphabet[node.sequence[0]]
            if node in nuc_seqs:
                node.sequence = nuc_seqs[node]
            if node.up is not None:
                node.__setattr__(attr+'_transitions', node.mutations)
                if node in nuc_muts:
                    node.mutations = nuc_muts[node]

        self.tt.use_mutation_length=tmp_use_mutation_length

    def get_attr_list(self, get_attr):
        states = []
        for node in self.tree.find_clades():
            if get_attr in node.attr:
                states.append(node.attr[get_attr])
        return states

    def add_translations(self):
        '''
        translate the nucleotide sequence into the proteins specified
        in self.proteins. these are expected to be SeqFeatures
        '''
        from Bio import Seq
        for node in self.tree.find_clades(order='preorder'):
            if not hasattr(node, "translations"):
                node.translations={}
                node.aa_mutations = {}
            if node.up is None:
                for prot in self.proteins:
                    node.translations[prot] = Seq.translate(str(self.proteins[prot].extract(Seq.Seq("".join(node.sequence)))).replace('-', 'N'))
                    node.aa_mutations[prot] = []
            else:
                for prot in self.proteins:
                    node.translations[prot] = Seq.translate(str(self.proteins[prot].extract(Seq.Seq("".join(node.sequence)))).replace('-', 'N'))
                    node.aa_mutations[prot] = [(a,pos,d) for pos, (a,d) in
                                    enumerate(zip(node.up.translations[prot],
                                                  node.translations[prot])) if a!=d]
        self.dump_attr.append('translations')


    def refine(self):
        '''
        add attributes for export, currently this is only muts and aa_muts
        '''
        self.tree.ladderize()
        for node in self.tree.find_clades():
            if node.up is not None:
                node.muts = ["".join(map(str, [a, pos+1, d])) for a,pos,d in node.mutations]
                node.aa_muts = {}
                if hasattr(node, 'translations'):
                    for prot in node.translations:
                        node.aa_muts[prot] = ["".join(map(str,[a,pos+1,d])) for a,pos,d in node.aa_mutations[prot]]
        for node in self.tree.find_clades(order="preorder"):
            if node.up is not None: #try:
                node.attr["div"] = node.up.attr["div"]+node.mutation_length
            else:
                node.attr["div"] = 0
        self.dump_attr.extend(['muts', 'aa_muts', 'aa_mutations', 'mutation_length', 'mutations'])


    def layout(self):
        """Add clade, xvalue, yvalue, mutation and trunk attributes to all nodes in tree"""
        clade = 0
        yvalue = self.tree.count_terminals()
        for node in self.tree.find_clades(order="preorder"):
            node.clade = clade
            clade += 1
            if node.up is not None: #try:
                node.xvalue = node.up.xvalue+node.mutation_length
                if self.is_timetree:
                    node.tvalue = node.numdate - self.tree.root.numdate
                else:
                    node.tvalue = 0
            else:
                node.xvalue = 0
                node.tvalue = 0
            if node.is_terminal():
                node.yvalue = yvalue
                yvalue -= 1
        for node in self.tree.get_nonterminals(order="postorder"):
            node.yvalue = np.mean([x.yvalue for x in node.clades])
        self.dump_attr.extend(['yvalue', 'xvalue', 'clade'])
        if self.is_timetree:
            self.dump_attr.extend(['tvalue'])


    def export(self, path = '', extra_attr = ['aa_muts', 'clade'], plain_export = 10):
        '''
        export the tree data structure along with the sequence information as
        json files for display in web browsers.
        parameters:
            path    -- path (incl prefix) to which the output files are written.
                       filenames themselves are standardized  to *tree.json and *sequences.json
            extra_attr -- attributes of tree nodes that are exported to json
            plain_export -- store sequences are plain strings instead of
                            differences to root if number of differences exceeds
                            len(seq)/plain_export
        '''
        from Bio import Seq
        from itertools import izip
        timetree_fname = path+'tree.json'
        sequence_fname = path+'sequences.json'
        tree_json = tree_to_json(self.tree.root, extra_attr=extra_attr)
        write_json(tree_json, timetree_fname, indent=None)

        # prepare a json with sequence information to export.
        # first step: add the sequence & translations of the root as string
        elems = {}
        elems['root'] = {}
        elems['root']['nuc'] = "".join(self.tree.root.sequence)
        for prot,seq in self.tree.root.translations.iteritems():
            elems['root'][prot] = seq

        # add sequence for every node in tree. code as difference to root
        # or as full strings.
        for node in self.tree.find_clades():
            if hasattr(node, "clade"):
                elems[node.clade] = {}
                # loop over proteins and nucleotide sequences
                for prot, seq in [('nuc', "".join(node.sequence))]+node.translations.items():
                    differences = {pos:state for pos, (state, ancstate) in
                                enumerate(izip(seq, elems['root'][prot]))
                                if state!=ancstate}
                    if plain_export*len(differences)<=len(seq):
                        elems[node.clade][prot] = differences
                    else:
                        elems[node.clade][prot] = seq

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
