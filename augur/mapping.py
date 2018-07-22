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

    prots = {}

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
                line = line.strip()
                dat = line.split('\t')
                POS = int(dat[posLoc])
                REF = dat[refLoc]
                ALT = dat[altLoc].split(',')
                GEN = dat[0] #'CHROM' or the gene name here
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
                    gen = GEN       #from CHROM, gene name

                    if gen not in prots.keys():
                        prots[gen] = {}
                        prots[gen]['sequences'] = {}
                        prots[gen]['positions'] = []
                        prots[gen]['reference'] = ''
                    if seq not in prots[gen]['sequences'].keys():
                        prots[gen]['sequences'][seq] = {}

                    #will never be insertion or deletion! because translation.
                    prots[gen]['sequences'][seq][pos] = alt
                    prots[gen]['positions'].append(pos)

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
        posN = np.array(prots[refSeq.name]['positions'])
        posN = np.unique(posN)
        prots[refSeq.name]['positions'] = list(np.sort(posN))

    return prots

def read_in_AA(drm_file):
    '''
    Reads in and stores position, alt AA, drug, and gene
    of mutations such as drug resistance mutations

    GENE	SITE	ALT_AA	DRUG
    gyrB	461	N	Fluoroquinolones
    gyrB	499	D	Fluoroquinolones
    rpoB	170	F	Rifampicin
    rpoB	359	A	Rifampicin

    '''
    import pandas as pd

    MUTs = defaultdict()

    df = pd.read_csv(drm_file, sep='\t')
    for mi, m in df.iterrows():
        pos = int(m.SITE)-1 #store in python numbering

        if m.GENE in MUTs:
            if pos not in MUTs[m.GENE]:
                MUTs[m.GENE][pos] = {}
                MUTs[m.GENE][pos]['base'] = [m.ALT_AA]
                MUTs[m.GENE][pos]['drug'] = m.DRUG
            else: 
                MUTs[m.GENE][pos]['base'].append(m.ALT_AA)
        else:
            MUTs[m.GENE] = {}
            MUTs[m.GENE][pos] = {}
            MUTs[m.GENE][pos]['base'] = [m.ALT_AA]
            MUTs[m.GENE][pos]['drug'] = m.DRUG

    site_info = MUTs

    return site_info


def read_in_nucs(drm_file):
    '''
    Reads in and stores position, alt base, drug, AA change (optional)
    of mutations such as drug-resistance mutations

    Format:

    GENOMIC_POSITION	ALT_BASE	AA	DRUG
    6505	T	D461N	Fluoroquinolones
    6505	C	D461N	Fluoroquinolones
    760314	T	V170F	Rifampicin
    760882	C	V359A	Rifampicin

    or:

    GENOMIC_POSITION	ALT_BASE	DRUG
    6505	T	Fluoroquinolones
    6505	C	Fluoroquinolones
    760314	T	Rifampicin
    760882	C	Rifampicin

    '''
    import pandas as pd

    MUTs = {}
    mutPositions = []

    df = pd.read_csv(drm_file, sep='\t')
    for mi, m in df.iterrows():
        pos = m.GENOMIC_POSITION-1 #put in python numbering
        mutPositions.append(pos)

        if pos in MUTs:
            MUTs[pos]['base'].append(m.ALT_BASE)
            if hasattr(m, "AA"):
                MUTs[pos]['AA'].append(m.AA)
        else:
            MUTs[pos] = {}
            MUTs[pos]['base'] = [m.ALT_BASE]
            MUTs[pos]['drug'] = m.DRUG
            if hasattr(m, "AA"):
                MUTs[pos]['AA'] = [m.AA]

    mutPositions = np.array(mutPositions)
    mutPositions = np.unique(mutPositions)
    mutPositions = np.sort(mutPositions)

    site_info = {'MUTs': MUTs,
                'mutPositions': mutPositions}

    return site_info

def find_aa(site_info, prots):
    '''
    Looks for DRM mutations which match in position and alt base in
    the translated protein dict

    '''
    seqMuts = {}

    for gene, val in site_info.items():
        for pos, info in val.items():
            if gene in prots and pos in prots[gene]['positions']:
                for seq, muts in prots[gene]['sequences'].items():
                    if pos in muts and muts[pos] in info['base']:
                        if seq not in seqMuts:
                            seqMuts[seq] = {}
                        refB = prots[gene]['reference'][pos]
                        seqMuts[seq][gene+": "+refB+str(pos+1)+muts[pos]] = info['drug']
                        #seqMuts[seq][refB+str(pos+1)+muts[pos]] = info['drug']  #if don't want to show genes
    return seqMuts


def find_sites(site_info, sequences, ref):
    '''
    Looks for DRM mutations which match in position and alt base in
    the sequence dict

    '''
    mutPositions = site_info['mutPositions']
    MUTs = site_info['MUTs']
    seqMuts = {}
    for key in mutPositions:
        for seq,v in sequences.items():
            try:
                base = sequences[seq][key]
                if base in MUTs[key]['base']:
                    if seq not in seqMuts:
                        seqMuts[seq] = {}
                    i=0
                    while base != MUTs[key]['base'][i]:
                        i+=1
                    if "AA" in MUTs[key]:
                        seqMuts[seq][MUTs[key]['AA'][i]] = MUTs[key]['drug']
                    else:
                        seqMuts[seq][ ref[key]+str(key+1)+MUTs[key]['base'][i] ] = MUTs[key]['drug']

            except KeyError:
                continue
    return seqMuts


def attach_muts(seqMuts, sequences, label, count):
    '''
    'Attaches' muts to nodes so format is correct to output as json
    '''
    mut_meta = defaultdict(lambda: {label: '0' })

    #Strings are used here because 0's dont work in auspice at moment
    #TODO change to ints when they do

    traitMuts = {}
    traitMuts[label] = ['0'] 
    for seq in sequences.keys():
        if seq in seqMuts:
            muts = 0
            for mut,drug in seqMuts[seq].items():
                drugs = drug.split(' ')
                for drug in drugs:
                    muts += 1
                    trDrug = drug
                    if trDrug in mut_meta[seq]:
                        mut_meta[seq][trDrug] = ",".join([mut_meta[seq][trDrug],mut])
                    else:
                        mut_meta[seq][trDrug] = mut

                    if trDrug in traitMuts:
                        if mut_meta[seq][trDrug] not in traitMuts[trDrug]:
                            traitMuts[trDrug].append(mut_meta[seq][trDrug])
                    else:
                        traitMuts[trDrug] = [ mut_meta[seq][trDrug] ]

            if count == "traits":
                numResist = str(len(set(mut_meta[seq].keys()))-1)
            else:
                numResist = str(muts)
            mut_meta[seq][label] = numResist
            if numResist not in traitMuts[label]:
                traitMuts[label].append(numResist)
        else:
            mut_meta[seq][label] = '0'

    return mut_meta, traitMuts

def run(args):
    '''
    This should be modified to work on Fasta-input files!!
    '''

    ## check file format and read in sequences
    is_vcf = False
    if any([args.alignment.lower().endswith(x) for x in ['.vcf', '.vcf.gz']]):
        if not args.vcf_reference:
            print("ERROR: a reference Fasta is required with VCF-format alignments")
            return -1
        is_vcf = True
        compress_seq = read_vcf(args.alignment, args.vcf_reference)
        sequences = compress_seq['sequences']
        ref = compress_seq['reference']
        aln = sequences
    else:
        # TO-DO fill in fasta-format processing
        aln = args.alignment

    #if give nucleotide positions:
    if args.nucleotide:
        site_info = read_in_nucs(args.nucleotide)

        #Find the MUTs and store them
        seqMuts = find_sites(site_info, sequences, ref)

    elif args.amino_acid:
        site_info = read_in_AA(args.amino_acid)

        prots = read_in_translate_vcf(args.translate_vcf, args.translate_ref)
        seqMuts = find_aa(site_info, prots)

    #Get json format; collect all unique options to generate colours for
    mut_meta, traitMuts = attach_muts(seqMuts, sequences, args.label, args.count)

    #write out json
    with open(args.output, 'w') as results:
        json.dump({"nodes":mut_meta}, results, indent=1)
