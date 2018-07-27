import numpy as np
from treetime.vcf_utils import read_vcf
from collections import defaultdict
import json

def read_in_translate_vcf(vcf_file, ref_file):
    """
    Reads in a vcf file where TRANSLATIONS have been stored and associated
    reference sequence fasta (to which the VCF file is mapped)
    This is the file output by "write_VCF_translation" below

    Very simple compared to the above as will never be insertion or deletion

    Returns a nested dict in the same format as is *input* in "write_VCF_translation" below,
    with a nested dict for each gene, which contains 'sequences', 'positions', and 'reference'
    """
    from Bio import SeqIO
    import numpy as np

    def mutation_struct():
        return {'sequences':defaultdict(dict), 'positions':[], 'reference':''}
    prots = defaultdict(mutation_struct)

    posLoc = 0
    refLoc = 0
    altLoc = 0
    sampLoc = 9

    #Use different openers depending on whether compressed
    opn = gzip.open if vcf_file.endswith(('.gz', '.GZ')) else open

    with opn(vcf_file, mode='rt') as f:
        for line in f:
            if line[0] != '#':
                #actual data
                dat = line.strip().split('\t')
                POS = int(dat[posLoc])
                REF = dat[refLoc]
                ALT = dat[altLoc].split(',')
                GENE = dat[0] #'CHROM' or the gene name here
                calls = np.array(dat[sampLoc:])

                #get samples that differ from Ref at this site
                recCalls = {}
                for sname, sa in zip(samps, calls):
                    if sa != '.':
                        recCalls[sname] = sa

                #store position and the altLoc
                for seq, gen in recCalls.items():
                    alt = str(ALT[int(gen[0])-1])   #get the index of the alternate
                    ref = REF
                    pos = POS-1     #VCF numbering starts from 1, but Reference seq numbering
                                    #will be from 0 because it's python!
                    gene = GENE       #from CHROM, gene name

                    #will never be insertion or deletion! because translation.
                    prots[gene]['sequences'][seq][pos] = alt
                    prots[gene]['positions'].append(pos)

            elif line[0] == '#' and line[1] == 'C':
                #header line, get all the information
                header = line.strip().split('\t')
                posLoc = header.index("POS")
                refLoc = header.index("REF")
                altLoc = header.index("ALT")
                sampLoc = header.index("FORMAT")+1
                samps = header[sampLoc:]
                nsamp = len(samps)

    for refSeq in SeqIO.parse(ref_file, format='fasta'):
        prots[refSeq.name]['reference'] = str(refSeq.seq)
        posN = np.unique(prots[refSeq.name]['positions'])
        prots[refSeq.name]['positions'] = list(posN)

    return prots


def read_in_features(drm_file):
    '''
    Reads in and stores position, alt base, feature, AA change (optional)
    of mutations such as drug-resistance mutations

    Format:

    SITE	ALT	DISPLAY_NAME	FEATURE
    6505	T	D461N	Fluoroquinolones
    6505	C	D461N	Fluoroquinolones
    760314	T	V170F	Rifampicin
    760882	C	V359A	Rifampicin

    or:

    SITE	ALT	FEATURE
    6505	T	Fluoroquinolones
    6505	C	Fluoroquinolones
    760314	T	Rifampicin
    760882	C	Rifampicin

    '''
    import pandas as pd

    MUTs = defaultdict(lambda:defaultdict(dict))

    mutPositions = defaultdict(list)

    df = pd.read_csv(drm_file, sep='\t' if drm_file.endswith('.tsv') else ',')
    for mi, m in df.iterrows():
        pos = m.SITE-1 #put in python numbering
        gene = m.GENE if hasattr(m, 'GENE') else 'nuc'

        mutPositions[gene].append(pos)
        MUTs[gene][pos][m.ALT] = {'feature':m.FEATURE.split()}
        if hasattr(m, "DISPLAY_NAME"):
            MUTs[gene][pos][m.ALT]['display_name'] = m.DISPLAY_NAME

    for gene in mutPositions:
        mutPositions[gene] = np.unique(mutPositions[gene])

    return MUTs


def annotate_strains_by_gene(annotations, features, sequences, gene='nuc'):
    for pos, info in features.items():   # annotated positions in gene
        # is the site mutated anywhere in the dataset?
        if pos in sequences['positions']:
            # loop over sequences in dataset
            for seq_name, muts in sequences['sequences'].items():
                # if position mutated to the right base/aa
                if pos in muts and muts[pos] in info:
                    der = muts[pos]
                    anc = sequences['reference'][pos]
                    feat = info[der]['feature']
                    if 'display_name' in info[der]:
                        label = info[der]['display_name']
                    elif gene=='nuc':
                        label = anc+str(pos+1)+der
                    else:
                        label = gene+": "+anc+str(pos+1)+der
                    annotations[seq_name][label] = feat


def annotate_strains(all_features, all_sequences):
    '''
    Looks for DRM mutations which match in position and alt base in
    the translated protein dict

    '''
    annotations = defaultdict(dict)
    for gene, features in all_features.items(): # loop over genes
        if gene in all_sequences:
            # this updates annotations in place!
            annotate_strains_by_gene(annotations, features, all_sequences[gene], gene=gene)

    return annotations


def attach_features(annotations, sequences, label, count):
    '''
    'Attaches' muts to nodes so format is correct to output as json
    '''

    #Strings are used here because 0's dont work in auspice at moment
    #TODO change to ints when they do
    seq_feature_dict = defaultdict(lambda: {label: '0' })

    for seq, anno in annotations.items():
        muts = 0
        for mut, features in anno.items():
            for feat in features:
                muts += 1
                if feat in seq_feature_dict[seq]:
                    seq_feature_dict[seq][feat] += ","+str(mut)
                else:
                    seq_feature_dict[seq][feat] = mut

        if count == "traits":
            numResist = str(len(set(seq_feature_dict[seq].keys()))-1)
        else:
            numResist = str(muts)

        seq_feature_dict[seq][label] = numResist

    return seq_feature_dict


def run(args):
    '''
    This should be modified to work on Fasta-input files!!
    '''

    ## check file format and read in sequences
    is_vcf = False
    if any([args.alignment.lower().endswith(x) for x in ['.vcf', '.vcf.gz']]):
        if not args.vcf_reference:
            print("ERROR: a reference Fasta is required with VCF-format alignments")
            return 1
        is_vcf = True
        if args.amino_acid:
            compress_seq = read_in_translate_vcf(args.alignment, args.vcf_reference)
        else:
            compress_seq = {"nuc": read_vcf(args.alignment, args.vcf_reference)}
    else:
        # TO-DO fill in fasta-format processing
        aln = args.alignment

    features = read_in_features(args.features)
    annotations = annotate_strains(features, compress_seq)
    #convert the annotations into string label that auspice can display
    seq_features = attach_features(annotations, compress_seq, args.label, args.count)

    #write out json
    with open(args.output, 'w') as results:
        json.dump({"nodes":seq_features}, results, indent=1)
