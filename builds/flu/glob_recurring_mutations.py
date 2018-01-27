from __future__ import print_function
import os, sys, glob, json

lineages = ['H3N2', 'H1N1pdm', 'Vic', 'Yam']
res = ['2y' ,'3y', '6y']
passage = ['cell', 'egg']
assay = ['hi', 'fra']

def return_struct():
    return {k:{k2:{k3:[] for k3 in assay} for k2 in passage} for k in res}

recurring_positions = {'H3N2':return_struct(), 'H1N1pdm':return_struct(),
                    'Vic':return_struct(), 'Yam':return_struct()}

recurring_mutations = {'H3N2':return_struct(), 'H1N1pdm':return_struct(),
                    'Vic':return_struct(), 'Yam':return_struct()}

r='2y'
for l in lineages:
    for p in passage:
        for a in assay:
            try:
                fname = 'processed/recurring_mutations/flu_%s_%s_%s_%s_ha_recurring_positions.txt'%(l.lower(), r, p, a)
                with open(fname, 'r') as ifile:
                    tmp =[]
                    for line in ifile:
                        if line[0]=='#' and len(line.strip()):
                            continue
                        tmp.append(line.strip().split())
                recurring_positions[l][r][p][a] = tmp

                fname = 'processed/recurring_mutations/flu_%s_%s_%s_%s_ha_recurring_mutations.txt'%(l.lower(), r, p, a)
                with open(fname, 'r') as ifile:
                    tmp =[]
                    for line in ifile:
                        if line[0]=='#' and len(line.strip()):
                            continue
                        tmp.append(line.strip().split())
                recurring_mutations[l][r][p][a] = tmp
            except:
                pass

with open('auspice/recurring_mutations.json', 'w') as ofile:
    json.dump(recurring_mutations, ofile)

with open('auspice/recurring_positions.json', 'w') as ofile:
    json.dump(recurring_positions, ofile)

