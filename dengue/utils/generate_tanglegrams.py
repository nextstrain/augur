import subprocess
import sys
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio import AlignIO
import re
from glob import glob
# import colored_traceback.always
import argparse


#########     Parse input, define functions for alignment splitting and tree building     ############

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--alignment', default=None, type=str, help="Full path of genome alignment file")
parser.add_argument('-r', '--reference', default=None, type=str, help="Full path of genbank format reference file to pull gene coordinates from")
parser.add_argument('-p', '--proteins', default=['C', 'M', 'E', 'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', '2K', 'NS4B', 'NS5'], nargs='*', help="Names of proteins to include (must match genbank annotations). E.g.: `-p C M E NS1`")
parser.add_argument('-o', '--outgroup', type=str, default = None, help="Taxa name to use as outgroup when rooting tree")
parser.add_argument('-f', '--fastafield', type=int, default = 1, help="Index of | delimited fasta header field where accession numbers live. E.g., `-f 1` for >whatever|accession|blah")
parser.add_argument('--raxmlpath', type=str, default = 'raxml', help="Path to raxml. If raxml is already in your .bashrc, leave default (just calls 'raxml')")
parser.add_argument('-run', nargs='*', default=['a', 't', 'p'],  choices=['a', 't', 'p'], help="Which operations to do: align (a), build trees (t), plot tanglegrams (p). E.g.: `--run a t` will split the alignment, build trees, and quit. Alternatively, `--run p` will find the raxml bestTree files in the cwd and plot a tanglechain.")

args = parser.parse_args()
alignment, reference = args.alignment, args.reference
alignmentobject = AlignIO.read(open(alignment, 'r'), 'fasta')
protein_list = args.proteins # [ protein1, protein2, ... ]
proteins = None #[ (protein1, FeatureLocation), (protein2, FeatureLocation), ...]

def load_reference(reference_file):
    reference_seq = SeqIO.read(reference_file, 'genbank')
    seen = []
    genome_annotation = reference_seq.features
    proteins = {}
    for f in reference_seq.features:
        if 'gene' in f.qualifiers and f.qualifiers['gene'][0] in protein_list:
            proteins[f.qualifiers['gene'][0]] = FeatureLocation(start=f.location.start, end=f.location.end, strand=1)
            seen.append(f.qualifiers['gene'][0])
    proteins = sorted(proteins.items(), key=lambda i: protein_list.index(i[0]))
    return proteins, reference_seq

def convert_coordinate(refseq, compareseq, coordinate):
    '''
    In: ungapped reference sequence, gapped (aligned) compare sequence, position to convert
    Out: coordinate of corresponding condition in the aligned sequence
    '''
    coordinate = coordinate - 1 # Adjust for python coordinates
    #check to make sure we have at least 100bp of downstream sequence, no more than 40bp of which are gaps, to match on
    assert len(refseq[coordinate:]) >= 100 and refseq[coordinate:coordinate+100].count('-')<40, 'ERROR: Not enough downstream context available for '+str(coordinate+1)
    #check to make sure the given coordinate doesn't correspond to a gap in the reference sequence
    assert refseq[coordinate] != '-', 'ERROR! Coordinate '+str(coordinate+1)+' is a gap in the reference sequence.'

    reference = refseq[coordinate:coordinate+100].replace('-', '')
    refpattern = '-*'.join(list(reference)) # Match in new sequence while ignoring gaps
    matchlist = re.findall(refpattern, compareseq, flags=re.IGNORECASE) # Check to make sure we get one and only one match
    assert len(matchlist) == 1, 'ERROR: found %d matches for coordinate %d'%(len(matchlist), coordinate+1)
    return re.search(refpattern, compareseq, flags=re.IGNORECASE).start()+1 #return the converted coordinate (adjusted back to genomic coordinates)

def split_alignment(alignment_file, proteins):
    '''
    in: path to fasta-format alignment file, [(gene, (start,end)), ...]
    out: separate .phyx-format (relaxed phylip) alignment files for each gene, based on provided coordinates
    NB: removes any sequences without > 70% non-gap sites from each segment alignment
    '''
    align = AlignIO.read(open(alignment_file, 'r'), 'fasta')
    ofile_stem = alignment_file.split('/')[-1].split('.')[0] # name output like alignmentfilename_protein.phyx
    ofile_list = []
    for protein, (start, end) in proteins:
        ofile_name = ofile_stem+'_%s.phyx'%(protein)
        ofile_list.append(ofile_name)
        start, end = start - 1, end - 1 # adjust for pythonic coordinates
        align_segment = align[:, start:end+1] #[allrows, startcolumn:endcolumn] endcolumn += 1 for inclusive slicing
        filtered_align_segment = AlignIO.MultipleSeqAlignment([])
        for seq in align_segment:
            seq.seq.data = str(seq.seq).replace('n', '-').replace('N', '-')
            if float(str(seq.seq).count('-')) / float(len(str(seq.seq))) <= 0.30: # Require at least 70% of sites are not gaps to include in segment alignment
                filtered_align_segment.append(seq)
        AlignIO.write(filtered_align_segment, ofile_name, 'phylip-relaxed')
    return ofile_list

def run_raxml(alignment_list):
    '''
    in: alignment files in .phyx format
    out: raxml ML trees (see below for parameters)
    NB: removes all raxml-generated files other than bestTree
    '''
    for a in alignment_list:
        protein = a.split('_')[-1].split('.')[0]
        # rapid hill climbing method; 2 threads; GTR model with -c categories of rate variation; specify seed for parsimony (ensure reproducibility); in/out fnames.
        try:
            subprocess.check_call("%s -f d -T 2 -m GTRCAT -c 25 -p 344312987 -o '%s' -s %s -n topology_%s"%(args.raxmlpath, outgroup, a, protein), shell=True)
        except:
            print 'WARNING: OS Error raised. Most likely, your outgroup %s was not found in the alignment for protein %s. Rerunning raxml with no outgroup; tanglechain for this segment may look artificially extra tangled. If it fails again, you have some other raxml config errorself.\n\n'%(outgroup, protein)
            subprocess.call('rm RAxML*%s'%protein, shell=True)
            subprocess.check_call("%s -f d -T 2 -O -m GTRCAT -c 25 -p 344312987 -s %s -n topology_%s"%(args.raxmlpath, a, protein), shell=True)
        subprocess.call(['rm RAxML_info*', 'rm RAxML_checkpoint*', 'rm RAxML_parsimonyTree*', 'rm RAxML_result*', 'rm RAxML_log*', 'rm *.reduced'], shell=True) # Cleanup


############ Plotting Utilities ######################

def euclidean((x1,y1),(x2,y2)):
    '''Given two (x,y) points as input, returns euclidean distance between them.'''
    return ((float(x2)-float(x1))**2+(float(y1)-float(y2))**2)**0.5

def sum_tip_distances(tree1,tree2):
    '''
    In: two baltic tree Objects
    Out: summed plotting distance for each tip that is present in both trees
    '''
    tree1_tips = { k.numName: k for k in tree1.Objects if k.branchType=='leaf' }
    tree2_tips = { k.numName: k for k in tree2.Objects if k.branchType=='leaf' }
    shared_tips = set(tree1_tips.keys()).intersection(set(tree2_tips.keys()))
    total_dist = 0.0

    for t in shared_tips:
        total_dist += euclidean( (tree1_tips[t].x, tree1_tips[t].y), (tree2_tips[t].x, tree2_tips[t].y) )

    return total_dist

def untangle(tree1, tree2):
    '''
    In: two baltic tree objects
    Does: Brute force method; iterates through every node of second tree (starting at root),
    rotates it iif reduces the total plotting distance between the shared tips.
    '''
    current_distance = sum_tip_distances(tree1, tree2)

    for n in sorted(tree2.nodes,key=lambda x: -x.height):
        if n.parent=='Root':
            continue
        n.rotate()
        tree2.drawTree()
        new_distance = sum_tip_distances(tree1, tree2)
        if new_distance <= current_distance:
            current_distance = new_distance
            continue
        else:
            n.rotate()
            tree2.drawTree()

#######         Run         #########################
if 'a' in args.run:
    proteins, reference_seq = load_reference(reference) #[ (protein1, FeatureLocation), (protein2, FeatureLocation), ...], SeqIO.SequenceObject
    reference_acc, reference_seq = reference_seq.id.split('.')[0], str(reference_seq.seq) # NB: Expects reference_seq.id like accession.whateverversionignored
    print 'Splitting alignment for each protein:', [p[0] for p in proteins], '\n\n'

    compare_seq = None
    for i in SeqIO.parse(alignment, 'fasta'): # Find the reference sequence in the alignment by matching accession numbers.
        if i.description.split('|')[args.fastafield].split('.')[0] == reference_acc:
            compare_seq = str(i.seq)
    assert compare_seq != None

    converted_proteins = [] # Convert coordinates from reference to aligned
    for gene, loc in proteins: #[ (protein1, (start, end)), (protein2, (start, end)), ...]
        start = convert_coordinate(reference_seq, compare_seq, int(loc.start))
        end = convert_coordinate(reference_seq, compare_seq, int(loc.end))
        converted_proteins.append((gene, (start, end)))

    alignments = split_alignment(alignment, converted_proteins)

if 't' in args.run:
    if 'a' not in args.run: # If we didn't generate segment alignments this run, look for them in the cwd as .phyx files
        alignments = [ t for t in glob('*.phyx') if t.split('_')[-1].split('.')[0] in protein_list ]
        print 'Building trees for %d alignments found in the cwd:\n'%len(alignments), alignments, '\n\n'
    else:
        print 'Building trees for alignments:\n', alignments, '\n\n'
    assert len(alignments) > 1

    if args.outgroup != None: # If no outgroup taxa provided, try to use the reference sequence.
        print 'Using specified outgroup to build trees.\n\n'
        outgroup = args.outgroup
    elif reference != None:
        print 'No outgroup provided; using provided reference sequence as outgroup.\n\n'
        outgroup = load_reference(reference)[1].name
    else:
        raise ValueError, 'Must provide outgroup or reference sequence name to root trees.\n\n'

    run_raxml(alignments) # Build trees.

if 'p' not in args.run:
    print 'Not plotting tanglechain per -run argument.\n\n'
    sys.exit()
else:
    pass

######## Prepare for plotting   #############
from matplotlib import pyplot as plt
import matplotlib as mpl
try:
    import baltic as bt
except Exception as e:
    print 'The baltic module is available from https://github.com/blab/baltic'
    raise e

print 'Plotting tanglechain in this order:\n', protein_list, '\n\n'

try:
    treefiles = sorted([ t for t in glob('*bestTree*') if t.split('_')[-1] in protein_list], key=lambda t: protein_list.index(t.split('_')[-1]))
    assert len(treefiles) == len(protein_list)
except ValueError as e: # For now, require that we can match all tree files to a protein
    print 'Oops! Plotting without building trees? Make sure your trees are named like `whatever_protein` and you passed a list of matching `protein` names to `-p` to set the order of trees.\n\n'
    raise e
except AssertionError as e:
    print 'ERROR: Missing tree files. Looked for trees for these proteins:\n', protein_list, '\n\n'
    raise e

trees = {}
for i, t in enumerate(treefiles):
    treestring, treeobject = open(t, 'r').readline().strip(), bt.tree()
    bt.make_tree(treestring, treeobject)
    treeobject.treeStats() ## initial traversal, checks for stats
    treeobject.sortBranches() ## traverses tree, sorts branches, draws tree (sets plotting coordinates)
    trees[i] = treeobject

for i in range(1,len(treefiles)):
    print 'Untangling tree number %d'%i
    untangle(trees[i-1], trees[i])

################
## Plot Genome Map
################
if proteins == None and reference != None: # If we didn't parse proteins earlier, but have the reference sequence, do so now.
    proteins, reference_seq = load_reference(reference)
    reference_seq = str(reference_seq.seq)

if proteins != None:                        # If possible, plot a to-scale genome map below the tanglechain.
    fig, (ax, genome) = plt.subplots(2, figsize=(25,15), gridspec_kw = {'height_ratios':[6, 1]})
    proteins = [ (p, (int(loc.start), int(loc.end))) for (p, loc) in proteins ]
    genomeL = len(reference_seq)

    def rescale(x):
        return (float(x)/genomeL)

    for gene, i in proteins:
        length = rescale(i[1]-i[0])
        genome.text(rescale(i[0])+0.5*length, 0.1,'%s'%(gene),va='center',ha='center',size=28,zorder=11)

        c = 'lightgray'
        genome.arrow(rescale(i[0]), 0.1, length, 0.0, alpha=0.6,head_width=0.15, width=0.25,head_length=0.0,length_includes_head=True,facecolor=c)

    genome.set_ylim(0,1) ## set y limits
    genome.set_xlim(0, 1)

    genome.spines['top'].set_visible(False)
    genome.spines['right'].set_visible(False)
    genome.spines['left'].set_visible(False)
    genome.spines['bottom'].set_visible(False)

    genome.tick_params(axis='x',labelsize=0,size=0)
    genome.tick_params(axis='y',labelsize=0,size=0)

else:                                   # If we have no annotations, skip the genome plot and just plot trees
    fig,ax = plt.subplots(figsize=(25,10))

######################
## Plot Trees
#####################

tree_names=protein_list ## define order in which dict will be accessed
tip_positions={x:{} for x in tree_names} ## remember the position of each tip in each tree

traitName='None' ## choose a trait to colour branches by
cmap=mpl.cm.viridis
cumulative_displace=0 ## this tracks the "current" x position, so trees are plotted one after another

branchWidth=2 ## increase branch width, since trees will be smaller
gapsize = 0.3*max([tree.treeHeight for tree in trees.values()])

for t, tr in enumerate(tree_names): ## iterate over trees
    cur_tree=trees[t] ## fetch tree object
    for k in cur_tree.Objects: ## iterate over branches
        if isinstance(k,bt.leaf): ## only interested in leaves
            tip_positions[tr][k.numName]=(k.height,k.y) ## remember tree, tip's position

for t,tr in enumerate(tree_names): ## iterate over trees
    cur_tree=trees[t] ## fetch tree object

    for k in cur_tree.Objects: ## iterate over branches
        x=k.x ## or from x position determined earlier
        y=k.y ## get y position from .drawTree that was run earlier, but could be anything else

        xp=k.parent.x ## get x position of current object's parent
        if x==None: ## matplotlib won't plot Nones, like root
            x=0.0
        if xp==None:
            xp=x

        x+=cumulative_displace ## adjust branch position by displacement, which depends on the position of tree in the overall plot
        xp+=cumulative_displace ## same for branch's parent
        c=cmap(k.height/cur_tree.treeHeight) ## or be a function of something else
        if isinstance(k,bt.leaf): ## if leaf...
            s=30 ## tip size can be fixed
            try:
                pos_in_first_tree=tip_positions[tree_names[0]][k.numName][1] ## fetch y coordinate of same tip in the first tree
                frac_pos=pos_in_first_tree/float(len(cur_tree.Objects))*2.0 ## normalize coordinate to be within interval [0.0,1.0]

                ax.scatter(x,y,s=s,facecolor=cmap(frac_pos),edgecolor='none',zorder=11) ## plot circle for every tip
                ax.scatter(x,y,s=s+0.8*s,facecolor='k',edgecolor='none',zorder=10) ## plot black circle underneath

                if t!=len(tree_names)-1: ## as long as we're not at the last tree - connect tips with coloured lines
                    next_x,next_y=tip_positions[tree_names[t+1]][k.numName] ## fetch coordinates of same tip in next tree
                    next_x+=cumulative_displace+cur_tree.treeHeight+gapsize ## adjust x coordinate by current displacement and future displacement

                    ax.plot([x,next_x],[y,next_y],lw=1,ls='-',color=cmap(frac_pos),zorder=0) ## connect current tip with same tip in the next tree
            except:
                continue

        elif isinstance(k,bt.node): ## if node...
            ax.plot([x,x],[k.children[-1].y,k.children[0].y],lw=branchWidth,color='k',ls='-',zorder=9) ## plot vertical bar

        ax.plot([xp,x],[y,y],lw=branchWidth,color='k',ls='-',zorder=9) ## always plot branch

    cumulative_displace+=cur_tree.treeHeight+gapsize ## increment displacement by the height of the tree

ax.set_ylim(0.0,len([x for x in cur_tree.Objects if isinstance(x,bt.leaf)])+2.5*gapsize) ## set y limits
ax.set_xlim(-2.5*gapsize,cumulative_displace+2.5*gapsize)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

ax.tick_params(axis='x',labelsize=0,size=0)
ax.tick_params(axis='y',labelsize=0,size=0)

plt.savefig('tanglechain.png')
