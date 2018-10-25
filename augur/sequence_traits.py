"""
Annotate sequences based on amino-acid or nucleotide signatures.
"""

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

    Parameters
    ----------
    vcf_file : str
        name of the vcf file to be read, can be gzipped

    ref_file : str
        name of the fasta file with the reference sequence

    Returns
    -------
    dict
        dictionary of dictionaries with mutations of each strain for each sequence
        relative to the reference

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
    Reads in and stores position, alt base/AA, feature, gene,
    and 'display name' (optional) of mutations such 
    as drug-resistance mutations

    Format to map by both nucleotide and AA sites:
    ----------------------------------------------
    GENE	SITE	ALT	DISPLAY_NAME	FEATURE
    gyrB	461	N		Fluoroquinolones
    nuc	1472358	T	rrs: C513T	Streptomycin
    nuc	1673425	T	fabG1: C-15T	Isoniazid Ethionamide
    ethA	175	T		Ethionamide

    Format to map by AA site:
    -------------------------
    GENE	SITE	ALT	FEATURE
    gyrB	461	N	Fluoroquinolones
    gyrB	499	D	Fluoroquinolones
    rpoB	170	F	Rifampicin
    rpoB	359	A	Rifampicin

    Format to map by nucleotide site:
    -----------------------------------
    SITE    ALT    DISPLAY_NAME    FEATURE
    6505    T    D461N    Fluoroquinolones
    6505    C    D461N    Fluoroquinolones
    760314    T    V170F    Rifampicin
    760882    C    V359A    Rifampicin

    or to map by nucleotide site and display mutations:
    ---------------------------------------------------
    SITE    ALT    FEATURE
    6505    T    Fluoroquinolones
    6505    C    Fluoroquinolones
    760314    T    Rifampicin
    760882    C    Rifampicin

    Parameters
    ----------
    drm_file : str
        file defining sequence features to be used for annotations

    Returns
    -------
    dict
        dict of dict with sequence features index by gene name, position, and character state

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
        if hasattr(m, "DISPLAY_NAME") and not m.isnull().DISPLAY_NAME:
            MUTs[gene][pos][m.ALT]['display_name'] = m.DISPLAY_NAME

    for gene in mutPositions:
        mutPositions[gene] = np.unique(mutPositions[gene])

    return MUTs


def annotate_strains_by_gene(annotations, features, sequences, gene='nuc'):
    """Sort through all potential features and link them up with mutations to
    produce an annotation

    Parameters
    ----------
    annotations : dict
        dictionary of sequence features as read in by `read_in_features`.
        This is modified in place
    features : dict
        dictionary of features in one gene
    sequences : dict
        sequences of that gene
    gene : str, optional
        name of the gene
    """
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
    #need to record those with no mutations so they can be given zero counts later
    for seq_name in sequences['sequences'].keys():
        if seq_name not in annotations:
            annotations[seq_name] = {}


def annotate_strains(all_features, all_sequences):
    '''
    Looks for DRM mutations which match in position and alt base in
    the translated protein dict

    Parameters
    ----------
    all_features : dict
        dict of all features in all genes, will be processed gene by gene
    all_sequences : dict
        sequence dict of all genes

    Returns
    -------
    dict
        annotations based on the strains for each feature

    '''
    annotations = defaultdict(dict)
    for gene, features in all_features.items(): # loop over genes
        if gene in all_sequences:
            # this updates annotations in place!
            annotate_strains_by_gene(annotations, features, all_sequences[gene], gene=gene)

    return annotations


def attach_features(annotations, label, count):
    '''
    'Attaches' features to nodes and lists the corresponding mutations
    as values, that is
    >>>
        {nodename:{"Resistance 1":"mut1,mut2", "Resistance 2":"mut1"}}

    Parameters
    ----------
    annotations : dict
        annotations fo stgrains as globed together by `annotate_strains`
    label : label
        label of the feature set as specified by as command line argument
    count : str
        if equal to traits, will count the number of distinct features that
        occur in the annotation, otherwise will count the total number of mutations

    Returns
    -------
    dict
        json/dict to export
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


def register_arguments(parser):
    parser.add_argument('--ancestral-sequences', type=str, help="nucleotide alignment to search for sequence traits in")
    parser.add_argument('--translations', type=str, help="AA alignment to search for sequence traits in (can include ancestral sequences)")
    parser.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the nucleotide VCF was mapped to')
    parser.add_argument('--vcf-translate-reference', type=str, help='fasta file of the sequence the translated VCF was mapped to')
    parser.add_argument('--features', type=str, help='file that specifies sites defining the features in a tab-delimited format "GENOMIC_POSITION ALT_BASE DRUG AA(optional)"')
    parser.add_argument('--count', type=str, choices=['traits','mutations'], default='traits', help='Whether to count traits (ex: # drugs resistant to) or mutations')
    parser.add_argument('--label', type=str, default="# Traits", help='How to label the counts (ex: Drug_Resistance)')
    parser.add_argument('--output', '-o', type=str, help='output json with sequence features')


def run(args):
    '''
    This should be modified to work on Fasta-input files!!
    '''
    print("This method may change in future! Please use 'augur sequence-traits -h' to check the latest options.")
    ## check file format and read in sequences
    is_vcf = False
    if ( (args.ancestral_sequences and any([args.ancestral_sequences.lower().endswith(x) for x in ['.vcf', '.vcf.gz']])) or
        (args.translations and any([args.translations.lower().endswith(x) for x in ['.vcf', '.vcf.gz']])) ):
        if ((args.ancestral_sequences and not args.vcf_reference) or
            (args.translations and not args.vcf_translate_reference)):
            print("ERROR: a reference Fasta is required with VCF-format alignments")
            return 1
        is_vcf = True
        compress_seq = defaultdict(dict)
        if args.translations:
            compress_seq = read_in_translate_vcf(args.translations, args.vcf_translate_reference)
        if args.ancestral_sequences:
            compress_seq["nuc"] = read_vcf(args.ancestral_sequences, args.vcf_reference)
    else:
        # TO-DO fill in fasta-format processing
        aln = args.alignment

    features = read_in_features(args.features)
    annotations = annotate_strains(features, compress_seq)
    #convert the annotations into string label that auspice can display
    seq_features = attach_features(annotations, args.label, args.count)

    #write out json
    with open(args.output, 'w') as results:
        json.dump({"nodes":seq_features}, results, indent=1, sort_keys = True)
