import numpy as np
import colorsys
from treetime.vcf_utils import read_vcf
from collections import defaultdict
from .utils import write_json
import shutil, json

def read_in_DRMs(drm_file):
    '''
    Reads in and stores position, ref & alt base, drug, AA change
    of drug resistance mutations (DRMs)
    This reads in the file-format below - don't know if standardized!

    GENOMIC_POSITION	ALT_BASE	SUBSTITUTION	DRUG
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
            DRMs[pos]['AA'].append(m.SUBSTITUTION)
        else:
            DRMs[pos] = {}
            DRMs[pos]['base'] = [m.ALT_BASE]
            DRMs[pos]['drug'] = m.DRUG
            DRMs[pos]['AA'] = [m.SUBSTITUTION]

    drmPositions = np.array(drmPositions)
    drmPositions = np.unique(drmPositions)
    drmPositions = np.sort(drmPositions)

    DRM_info = {'DRMs': DRMs,
                'drmPositions': drmPositions}

    return DRM_info

def drugTranslate(x):
    '''
    As above file format gives drugs as abbreviations, convert here.
    '''
    return {
        'RIF': 'Rifampicin',
        'FQ': 'Fluoroquinolones',
        'ETH': 'Ethionamide',
        'EMB': 'Ethambutol',
        'SM': 'Streptomycin',
        'PZA': 'Pyrazinamide',
        'CAP': 'Capreomycin',
        'INH': 'Isoniazid',
        'KAN': 'Kanamycin',
        'AK': 'Amikacin'
    }[x]


def find_drms(DRM_info, sequences):
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
                    seqDRM[seq][DRMs[key]['AA'][i]] = DRMs[key]['drug']

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
            tempList = []
            for mut,drug in seqDRM[seq].items():
                drugs = drug.split(';')
                for drug in drugs:
                    trDrug = drugTranslate(drug)
                    if trDrug in drm_meta[seq]:
                        drm_meta[seq][trDrug] = ",".join([drm_meta[seq][trDrug],mut])
                    else:
                        drm_meta[seq][trDrug] = mut

                    if trDrug in drugMuts:
                        if drm_meta[seq][trDrug] not in drugMuts[trDrug]:
                            drugMuts[trDrug].append(drm_meta[seq][trDrug])
                    else:
                        drugMuts[trDrug] = [ drm_meta[seq][trDrug] ]

                tempList.append(trDrug)
            numResist = str(len(tempList))
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
        #fill in fasta-format processing
        aln = args.alignment

    DRM_info = read_in_DRMs(args.drms)

    #Find the DRMs and store them
    seqDRM = find_drms(DRM_info, sequences)

    #Get json format; collect all unique options to generate colours for
    drm_meta, drugMuts = attach_drms(seqDRM, sequences)

    #write out json
    with open(args.output, 'w') as results:
        json.dump({"nodes":drm_meta}, results, indent=1)
