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
    Reads in and stores position, ref & alt base, drug, AA change
    of drug resistance mutations (DRMs)
    This reads in the file-format below - don't know if standardized!

    GENOMIC_POSITION	ALT_BASE	AA	DRUG
    6505	T	D461N	FQ
    6505	C	D461N	FQ
    760314	T	V170F	RIF
    760882	C	V359A	RIF

    '''
    import pandas as pd

    DRMs = defaultdict()

    df = pd.read_csv(drm_file, sep='\t')
    for mi, m in df.iterrows():
        pos = int(m.SITE)

        if m.GENE in DRMs:
            if pos not in DRMs[m.GENE]:
                DRMs[m.GENE][pos] = {}
                DRMs[m.GENE][pos]['base'] = [m.ALT_AA]
                DRMs[m.GENE][pos]['drug'] = m.DRUG
            else: 
                DRMs[m.GENE][pos]['base'].append(m.ALT_AA)
            
        #if (m.GENE,pos) in DRMs:
        #    DRMs[(m.GENE,pos)]['base'].append(m.ALT_AA)
        else:
            DRMs[m.GENE] = {}
            DRMs[m.GENE][pos] = {}
            DRMs[m.GENE][pos]['base'] = [m.ALT_AA]
            DRMs[m.GENE][pos]['drug'] = m.DRUG

    DRM_info = DRMs

    return DRM_info


def read_in_DRMs(drm_file):
    '''
    Reads in and stores position, ref & alt base, drug, AA change
    of drug resistance mutations (DRMs)
    This reads in the file-format below - don't know if standardized!

    GENOMIC_POSITION	ALT_BASE	AA	DRUG
    6505	T	D461N	FQ
    6505	C	D461N	FQ
    760314	T	V170F	RIF
    760882	C	V359A	RIF

    '''
    import pandas as pd

    DRMs = {}
    drmPositions = []

    df = pd.read_csv(drm_file, sep='\t')
    for mi, m in df.iterrows():
        pos = m.GENOMIC_POSITION-1 #put in python numbering
        drmPositions.append(pos)

        if pos in DRMs:
            DRMs[pos]['base'].append(m.ALT_BASE)
            if hasattr(m, "AA"):
                DRMs[pos]['AA'].append(m.AA)
        else:
            DRMs[pos] = {}
            DRMs[pos]['base'] = [m.ALT_BASE]
            DRMs[pos]['drug'] = m.DRUG
            if hasattr(m, "AA"):
                DRMs[pos]['AA'] = [m.AA]

    drmPositions = np.array(drmPositions)
    drmPositions = np.unique(drmPositions)
    drmPositions = np.sort(drmPositions)

    DRM_info = {'DRMs': DRMs,
                'drmPositions': drmPositions}

    return DRM_info

def find_aa_drm(DRM_info, prots):
    seqDRM = {}

    for gene, val in DRM_info.items():
        for pos, info in val.items():
            if gene in prots and pos in prots[gene]['positions']:
                #import ipdb; ipdb.set_trace()
                for seq, muts in prots[gene]['sequences'].items():
                    if pos-1 in muts and muts[pos-1] in info['base']:
                        if seq not in seqDRM:
                            seqDRM[seq] = {}
                        refB = prots[gene]['reference'][pos-1]
                        seqDRM[seq][refB+str(pos)+muts[pos-1]] = info['drug']

    #import ipdb; ipdb.set_trace()
    return seqDRM
    #why isn't finding katG 315 INH

def find_drms(DRM_info, sequences, ref):
    '''
    Looks for DRM mutations which match in position and alt base in
    the sequence dict

    '''
    drmPositions = DRM_info['drmPositions']
    DRMs = DRM_info['DRMs']
    seqDRM = {}
    for key in drmPositions:
        for seq,v in sequences.items():
            try:
                base = sequences[seq][key]
                if base in DRMs[key]['base']:
                    if seq not in seqDRM:
                        seqDRM[seq] = {}
                    i=0
                    while base != DRMs[key]['base'][i]:
                        i+=1
                    if "AA" in DRMs[key]:
                        seqDRM[seq][DRMs[key]['AA'][i]] = DRMs[key]['drug']
                    else:
                        seqDRM[seq][ ref[key]+str(key+1)+DRMs[key]['base'][i] ] = DRMs[key]['drug']

            except KeyError:
                continue
    return seqDRM


def attach_drms(seqDRM, sequences):
    '''
    'Attaches' DRMs to nodes so format is correct to output as json
    Also collects up all DRM options that colours will need to be generated for
    '''
    drm_meta = defaultdict(lambda: {"Drug_Resistance": '0' })

    drugMuts = {}
    drugMuts["Drug_Resistance"] = ['0'] #strings because makes it easier to work elsewhere
    for seq in sequences.keys():
        if seq in seqDRM:
            for mut,drug in seqDRM[seq].items():
                drugs = drug.split(' ')
                for drug in drugs:
                    trDrug = drug
                    if trDrug in drm_meta[seq]:
                        drm_meta[seq][trDrug] = ",".join([drm_meta[seq][trDrug],mut])
                    else:
                        drm_meta[seq][trDrug] = mut

                    if trDrug in drugMuts:
                        if drm_meta[seq][trDrug] not in drugMuts[trDrug]:
                            drugMuts[trDrug].append(drm_meta[seq][trDrug])
                    else:
                        drugMuts[trDrug] = [ drm_meta[seq][trDrug] ]

            numResist = str(len(set(drm_meta[seq].keys()))-1)
            drm_meta[seq]["Drug_Resistance"] = numResist
            if numResist not in drugMuts["Drug_Resistance"]:
                drugMuts["Drug_Resistance"].append(numResist)
        else:
            drm_meta[seq]["Drug_Resistance"] = '0'

    return drm_meta, drugMuts

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
        DRM_info = read_in_DRMs(args.nucleotide)

        #Find the DRMs and store them
        seqDRM = find_drms(DRM_info, sequences, ref)

    elif args.amino_acid:
        DRM_info = read_in_AA(args.amino_acid)

        prots = read_in_translate_vcf(args.translate_vcf, args.translate_ref)
        seqDRM = find_aa_drm(DRM_info, prots)
        import ipdb; ipdb.set_trace()

    #Get json format; collect all unique options to generate colours for
    drm_meta, drugMuts = attach_drms(seqDRM, sequences)

    #write out json
    with open(args.output, 'w') as results:
        json.dump({"nodes":drm_meta}, results, indent=1)
