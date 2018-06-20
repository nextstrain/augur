import numpy as np
import colorsys
from treetime.vcf_utils import read_vcf
from collections import defaultdict
from .utils import write_json
import shutil, json

def read_in_DRMs(drm_file):
    '''
    Reads in and stores position, ref & alt base, drug, AA change, and gene
    of drug resistance mutations (DRMs)
    This reads in the file-format below - don't know if standardized!

    KEY	GENOMIC_POSITION	REF_BASE	ALT_BASE	LOCUS	GENE	SUBSTITUTION	DRUG	SOURCE	COMMENT
    1	6505	G	T	Rv0005	gyrB	D461N	FQ	Malik_PLOS_ONE_2012
    2	6505	G	C	Rv0005	gyrB	D461N	FQ	Malik_PLOS_ONE_2012
    3	6618	G	A	Rv0005	gyrB	N499D	FQ	Malik_PLOS_ONE_2012
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
            DRMs[pos]['gene'] = m.GENE

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

def get_N_HexCol(n=5):
    '''
    Auto-generate colours for DRMs to prevent users having to do it

    modified from jhrf's answer on stackoverflow: (modified to py3 from ceprio's answer)
    https://stackoverflow.com/questions/876853/generating-color-ranges-in-python
    '''
    import colorsys

    HSV_tuples = [(x*1.0/n, 0.7, 0.7) for x in range(n)]
    hex_out = []
    for rgb in HSV_tuples:
        rgb = map(lambda x: int(x*255),colorsys.hsv_to_rgb(*rgb))
        hex_out.append('#%02x%02x%02x' % tuple(rgb))
        #hex_out.append("#"+"".join(map(lambda x: chr(x).encode('hex'),rgb)))
    return hex_out


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
    '''
    alignment <- vcf alignment
    vcf_reference <- vcf ref
    drm_file <- drm file
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

    DRM_info = read_in_DRMs(args.drm_file)

    #Find the DRMs and store them
    seqDRM = find_drms(DRM_info, sequences)

    #Get json format; collect all unique options to generate colours for
    drm_meta, drugMuts = attach_drms(seqDRM, sequences)

    #write out json
    with open(args.output, 'w') as results:
        json.dump({"nodes":drm_meta}, results, indent=1)

    #Generate colours
    newColAdds = []
    drugMuts["Drug_Resistance"].sort()
    for drug,muts in drugMuts.items():
        cols = get_N_HexCol(len(muts))
        drugCol = [ "\t".join([drug, muts[i], cols[i]]) for i in range(len(muts)) ]
        drugCol.append("")
        newColAdds = newColAdds + drugCol

    #copy and add to existing color file, if supplied:
    write_status = 'w'
    if args.color_defs:
        shutil.copyfile(args.color_defs, args.col_output)
        write_status = 'a'
        newColAdds = ['',''] + newColAdds #skips a line, looks nice

    with open(args.col_output, write_status) as the_file:
        the_file.write("\n".join(newColAdds))










