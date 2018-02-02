from __future__ import print_function
import os, sys, glob, json

lineages = ['H3N2', 'H1N1pdm', 'Vic', 'Yam']
res = ['2y' ,'3y', '6y']
region = ['NA', 'AS', 'EU']


def return_struct():
    return {k:{k2:[] for k2 in region} for k in res}

rising_mutations = {'H3N2':return_struct(), 'H1N1pdm':return_struct(),
                    'Vic':return_struct(), 'Yam':return_struct()}



for l in lineages:
    for reg in region:
        fname = 'processed/rising_mutations/flu_%s_2y_ha_%s_rising_mutations.txt'%(l.lower(), reg)
        with open(fname, 'r') as ifile:
            tmp =[]
            for line in ifile:
                if line[0]=='#' and len(line.strip()):
                    continue
                tmp.append(line.strip().split())
        rising_mutations[l][r][reg] = tmp

with open('auspice/rising_mutations.json', 'w') as ofile:
    json.dump(rising_mutations, ofile)

