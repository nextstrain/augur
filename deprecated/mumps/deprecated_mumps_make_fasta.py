from __future__ import print_function
import os, sys, re
from Bio import SeqIO
from pdb import set_trace
from pprint import pprint

crazies = {
    'X63707_Jeryl_Lynn_Vaccine_Strain/01.63': ("Jeryl_Lynn_Vaccine_Strain_X63707","Philadelphia","USA","01","63"),
    'KX853532_Itatiba.SP.BRA/43.15': ("KX853532","Itatiba","BRA","15","43")
}

def fix_names(name):
    name = name.strip('"').strip("'")
    p = re.compile('([0-9a-zA-Z]+)_([a-zA-Z]+)\.([a-zA-Z]+)/(\d+)\.(\d+)')
    try:
        g = p.search(name).groups()
    except AttributeError:
        g = crazies[name]
    year_prefix = 19
    if int(g[4])<18:
        year_prefix = 20
    header = "{}|{}|{}|{}{}-{}".format(g[0],g[1],g[2],year_prefix,g[4],g[3])
    pprint(name + "->" + header)
    return(">" + header)

def fix_seqs(seq):
    return seq.upper().replace('U','T')

if __name__=="__main__":
    # records = list(SeqIO.parse("test_input/284_SH_Jenn_renamed_for_BEAST.nex", "nexus"))
    with open("test_input/284_SH_Jenn_renamed_for_BEAST.nex", "r") as fh:
        raw = fh.readlines()
    aaa = next(x for x in enumerate(raw) if "matrix" in x[1])
    raw = raw[aaa[0]+1:len(raw)]
    raw = [x.strip().split("\t") for x in raw if len(x)>10]
    data = {fix_names(a[0]):fix_seqs(a[1]) for a in raw}
    with open("test_input/284_SH_Jenn.mfa", 'w') as fh:
        for k, v in data.iteritems():
            fh.write(k+"\n"+v+"\n")
