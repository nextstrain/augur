from __future__ import division, print_function
import os, time, sys
sys.path.insert(0,'../') # path.abspath(path.join(__file__ ,".."))
from io_util import make_dir, remove_dir, tree_to_json, write_json, myopen
from sequences import sequence_set
import numpy as np
import math
from subprocess import CalledProcessError, check_call, STDOUT
from pprint import pprint
from pdb import set_trace
try:
    from treetime_augur import TreeTime, utils
except ImportError:
    print("Couldn't import treetime_augur. Here's the searched paths:")
    pprint(sys.path)

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
    def __init__(self, aln, proteins=None, verbose=2, logger=None, **kwarks):
        super(tree, self).__init__()
        self.aln = aln
        self.nthreads = 2
        self.sequence_lookup = {seq.id:seq for seq in aln}
        self.nuc = kwarks['nuc'] if 'nuc' in kwarks else True
        self.dump_attr = []
        self.verbose = verbose
        if proteins!=None:
            self.proteins = proteins
        else:
            self.proteins={}
        if 'run_dir' not in kwarks:
            import random
            self.run_dir = '_'.join(['temp', time.strftime('%Y%m%d-%H%M%S',time.gmtime()), str(random.randint(0,1000000))])
        else:
            self.run_dir = kwarks['run_dir']
        if logger is None:
            def f(x,y):
                if y>self.verbose: print(x)
            self.logger = f
        else:
            self.logger=logger

    def getDateMin(self):
        return self.tree.root.date

    def getDateMax(self):
        dateMax = self.tree.root.date
        for node in self.tree.find_clades():
            if node.date > dateMax:
                dateMax = node.date
        return dateMax

    def dump(self, treefile, nodefile):
        from Bio import Phylo
        Phylo.write(self.tree, treefile, 'newick')
        node_props = {}
        for node in self.tree.find_clades():
            node_props[node.name] = {attr:node.__getattribute__(attr) for attr in self.dump_attr if hasattr(node, attr)}

        with myopen(nodefile, 'w') as nfile:
            from cPickle import dump
            dump(node_props, nfile)


    def build(self, root='midpoint', raxml=True, raxml_time_limit=1, raxml_bin='raxml', debug=False, num_distinct_starting_trees=1):
        from Bio import Phylo, AlignIO
        import subprocess, glob, shutil
        make_dir(self.run_dir)
        os.chdir(self.run_dir)
        for seq in self.aln: seq.name=seq.id
        AlignIO.write(self.aln, 'temp.fasta', 'fasta')
        out_fname = "tree_infer.newick"
        if raxml:
            self.logger("modified RAxML script - no branch length optimisation or time limit", 1)
            AlignIO.write(self.aln,"temp.phyx", "phylip-relaxed")
            if num_distinct_starting_trees == 1:
                cmd = raxml_bin + " -f d -T " + str(self.nthreads) + " -m GTRCAT -c 25 -p 235813 -n tre -s temp.phyx"
            else:
                self.logger("RAxML running with {} starting trees (longer but better...)".format(num_distinct_starting_trees), 1)
                cmd = raxml_bin + " -f d -T " + str(self.nthreads) + " -N " + str(num_distinct_starting_trees) + " -m GTRCAT -c 25 -p 235813 -n tre -s temp.phyx"
            fh = open("raxml.log", 'w')
            try:
                check_call(cmd, stdout=fh, stderr=STDOUT, shell=True)
                self.logger("RAXML COMPLETED.", 1)
            except CalledProcessError:
                self.logger("RAXML TREE FAILED - check {}/raxml.log".format(self.run_dir), 1)
                sys.exit(2)
            shutil.copy('RAxML_bestTree.tre', out_fname)

            # if raxml_time_limit>0:
            #     tmp_tree = Phylo.read('initial_tree.newick','newick')
            #     resolve_iter = 0
            #     resolve_polytomies(tmp_tree)
            #     while (not tmp_tree.is_bifurcating()) and (resolve_iter<10):
            #         resolve_iter+=1
            #         resolve_polytomies(tmp_tree)
            #     Phylo.write(tmp_tree,'initial_tree.newick', 'newick')
            #     AlignIO.write(self.aln,"temp.phyx", "phylip-relaxed")
            #     self.logger( "RAxML tree optimization with time limit %d hours"%raxml_time_limit,2)
            #     # using exec to be able to kill process
            #     end_time = time.time() + int(raxml_time_limit*3600)
            #     process = subprocess.Popen("exec " + raxml_bin + " -f d -T " + str(self.nthreads) +
            #                                 " -j -s temp.phyx -n topology -c 25 -m GTRCAT -p 344312987"
            #                                 " -t initial_tree.newick > raxml1_stdout", shell=True)
            #     while (time.time() < end_time):
            #         if os.path.isfile('RAxML_result.topology'):
            #             break
            #         time.sleep(10)
            #     process.terminate()
            #
            #     checkpoint_files = glob.glob("RAxML_checkpoint*")
            #     if os.path.isfile('RAxML_result.topology'):
            #         checkpoint_files.append('RAxML_result.topology')
            #     if len(checkpoint_files) > 0:
            #         last_tree_file = checkpoint_files[-1]
            #         shutil.copy(last_tree_file, 'raxml_tree.newick')
            #     else:
            #         shutil.copy("initial_tree.newick", 'raxml_tree.newick')
            # else:
            #     shutil.copy("initial_tree.newick", 'raxml_tree.newick')
            #
            # try:
            #     self.logger("RAxML branch length optimization", 2)
            #     os.system(raxml_bin + " -f e -T " + str(self.nthreads)
            #               + " -s temp.phyx -n branches -c 25 -m GTRGAMMA -p 344312987 -t raxml_tree.newick  > raxml2_stdout")
            #     shutil.copy('RAxML_result.branches', out_fname)
            # except:
            #     self.logger("RAxML branch length optimization failed", 1)
            #     shutil.copy('raxml_tree.newick', out_fname)
        else:
            tree_cmd = ["fasttree"]
            if self.nuc: tree_cmd.append("-nt")
            tree_cmd.extend(["temp.fasta","1>","initial_tree.newick", "2>", "fasttree_stderr"])
            os.system(" ".join(tree_cmd))
            shutil.copy('initial_tree.newick', out_fname)
        self.tt_from_file(out_fname, root)
        os.chdir('..')
        if not debug:
            remove_dir(self.run_dir)


    def tt_from_file(self, infile, root='best', nodefile=None):
        self.is_timetree=False
        self.logger('Reading tree from file '+infile,2)
        dates  =   {seq.id:seq.attributes['num_date']
                    for seq in self.aln if 'date' in seq.attributes}
        self.tt = TreeTime(dates=dates, tree=infile, gtr='Jukes-Cantor',
                            aln = self.aln, verbose=self.verbose, fill_overhangs=True)
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
            self.logger('reading node properties from file: '+nodefile,2)
            with myopen(nodefile, 'r') as infile:
                from cPickle import load
                node_props = load(infile)
            for n in self.tree.find_clades():
                if n.name in node_props:
                    for attr in node_props[n.name]:
                        n.__setattr__(attr, node_props[n.name][attr])
                else:
                    self.logger("No node properties found for "+n.name,2)


    def ancestral(self, **kwarks):
        self.tt.optimize_seq_and_branch_len(infer_gtr=True, **kwarks)
        self.dump_attr.append('sequence')
        for node in self.tree.find_clades():
            if not hasattr(node,'attr'):
                node.attr = {}


    def timetree(self, Tc=0.01, infer_gtr=True, reroot='best', resolve_polytomies=True,
                 max_iter=2, confidence=False, use_marginal=False, **kwarks):
        self.logger('estimating time tree...',2)
        if confidence and use_marginal:
            marginal = 'assign'
        else:
            marginal = confidence
        self.tt.run(infer_gtr=infer_gtr, root=reroot, Tc=Tc, do_marginal=marginal,
                    resolve_polytomies=resolve_polytomies, max_iter=max_iter, **kwarks)
        self.logger('estimating time tree...done',3)
        self.dump_attr.extend(['numdate','date','sequence'])
        to_numdate = self.tt.date2dist.to_numdate
        for node in self.tree.find_clades():
            if hasattr(node,'attr'):
                node.attr['num_date'] = node.numdate
            else:
                node.attr = {'num_date':node.numdate}
            if confidence:
                node.attr["num_date_confidence"] = sorted(to_numdate(node.marginal_inverse_cdf([0.05,0.95])))

        self.is_timetree=True


    def geo_inference(self, attr, missing='?', root_state=None, report_likelihoods=False):
        '''
        infer a "mugration" model by pretending each region corresponds to a sequence
        state and repurposing the GTR inference and ancestral reconstruction
        '''
        from treetime_augur import GTR
        # Determine alphabet
        places = set()
        for node in self.tree.find_clades():
            if hasattr(node, 'attr'):
                if attr in node.attr and attr!=missing:
                    places.add(node.attr[attr])
        if root_state is not None:
            places.add(root_state)

        # construct GTR (flat for now). The missing DATA symbol is a '-' (ord('-')==45)
        places = sorted(places)
        nc = len(places)
        if nc>180:
            self.logger("geo_inference: can't have more than 180 places!",1)
            return
        elif nc==1:
            self.logger("geo_inference: only one place found -- setting every internal node to %s!"%places[0],1)
            for node in self.tree.find_clades():
                node.attr[attr] = places[0]
                node.__setattr__(attr+'_transitions',[])
            return
        elif nc==0:
            self.logger("geo_inference: list of places is empty!",1)
            return

        # store previously reconstructed sequences
        nuc_seqs = {}
        nuc_muts = {}
        nuc_seq_LH = None
        if hasattr(self.tt.tree,'sequence_LH'):
            nuc_seq_LH = self.tt.tree.sequence_LH
        for node in self.tree.find_clades():
            if hasattr(node, 'sequence'):
                nuc_seqs[node] = node.sequence
            if hasattr(node, 'mutations'):
                nuc_muts[node] = node.mutations
                node.__delattr__('mutations')


        alphabet = {chr(65+i):place for i,place in enumerate(places)}
        sequence_gtr = self.tt.gtr
        myGeoGTR = GTR.custom(pi = np.ones(nc, dtype=float)/nc, W=np.ones((nc,nc)),
                              alphabet = np.array(sorted(alphabet.keys())))
        missing_char = chr(65+nc)
        alphabet[missing_char]=missing
        myGeoGTR.profile_map[missing_char] = np.ones(nc)
        alphabet_rev = {v:k for k,v in alphabet.iteritems()}

        # set geo info to nodes as one letter sequence.
        for node in self.tree.get_terminals():
            if hasattr(node, 'attr'):
                if attr in node.attr:
                    node.sequence=np.array([alphabet_rev[node.attr[attr]]])
            else:
                node.sequence=np.array([missing_char])
        for node in self.tree.get_nonterminals():
            node.__delattr__('sequence')
        self.tt.make_reduced_alignment()
        if root_state is not None:
            self.tree.root.split(n=1, branch_length=0.0)
            extra_clade = self.tree.root.clades[-1]
            extra_clade.name = "dummy_root_node"
            extra_clade.up = self.tree.root
            extra_clade.sequence = np.array([alphabet_rev[root_state]])
        # set custom GTR model, run inference
        self.tt._gtr = myGeoGTR
        # import pdb; pdb.set_trace()
        tmp_use_mutation_length = self.tt.use_mutation_length
        self.tt.use_mutation_length=False
        self.tt.infer_ancestral_sequences(method='ml', infer_gtr=False,
            store_compressed=False, pc=5.0, marginal=True, normalized_rate=False)

        if root_state is not None:
            self.tree.prune(extra_clade)
        # restore the nucleotide sequence and mutations to maintain expected behavior
        self.tt.geogtr = self.tt.gtr
        self.tt.geogtr.alphabet_to_location = alphabet
        self.tt._gtr = sequence_gtr
        self.dump_attr.append(attr)
        if hasattr(self.tt.tree,'sequence_LH'):
            self.tt.tree.geo_LH = self.tt.tree.sequence_LH
            self.tt.tree.sequence_LH = nuc_seq_LH
        for node in self.tree.find_clades():
            node.attr[attr] = alphabet[node.sequence[0]]
            if node in nuc_seqs:
                node.sequence = nuc_seqs[node]
            if node.up is not None:
                node.__setattr__(attr+'_transitions', node.mutations)
                if node in nuc_muts:
                    node.mutations = nuc_muts[node]
            # save marginal likelihoods if desired
            if report_likelihoods:
                node.attr[attr + "_entropy"] = sum([v * math.log(v+1E-20) for v in node.marginal_profile[0]]) * -1 / math.log(len(node.marginal_profile[0]))
                # javascript: vals.map((v) => v * Math.log(v + 1E-10)).reduce((a, b) => a + b, 0) * -1 / Math.log(vals.length);
                self.dump_attr.append(attr + "_entropy")
                marginal = [(alphabet[self.tt.geogtr.alphabet[i]], node.marginal_profile[0][i]) for i in range(0, len(self.tt.geogtr.alphabet))]
                marginal.sort(key=lambda x: x[1], reverse=True) # sort on likelihoods
                marginal = [(a, b) for a, b in marginal if b > 0.01][:4] #only take stuff over 1% and the top 4 elements
                self.dump_attr.append(attr + "_likelihoods")
                node.attr[attr + "_likelihoods"] = {a:b for a,b in marginal}
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


    def export(self, path = '', extra_attr = ['aa_muts', 'clade'], plain_export = 10, indent=None):
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
        timetree_fname = path+'_tree.json'
        sequence_fname = path+'_sequences.json'
        tree_json = tree_to_json(self.tree.root, extra_attr=extra_attr)
        write_json(tree_json, timetree_fname, indent=indent)

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

        write_json(elems, sequence_fname, indent=indent)


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
