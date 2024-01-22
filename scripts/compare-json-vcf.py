import argparse
from Bio import SeqIO
from treetime.vcf_utils import read_vcf
import sys, json
from augur.utils import read_tree
from collections import defaultdict

if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "Compare mutations in VCF vs JSON")
    parser.add_argument('--json', required=True, type=str, help="node-data JSON from augur ancestral")
    parser.add_argument('--vcf', required=True, type=str, help="VCF output from augur ancestral")
    parser.add_argument('--ref', required=True, type=str, help="reference fasta (supplied to augur ancestral)")
    parser.add_argument('--tree', required=True, type=str, help="newick tree (supplied to augur ancestral)")
    args = parser.parse_args()

    vcf = read_vcf(args.vcf, args.ref)
    ref = str(SeqIO.read(args.ref, 'fasta').seq).upper()
    with open(args.json) as fh:
        ndj = json.load(fh)
    T = read_tree(args.tree)
    mask = ndj['mask']

    # import pdb; pdb.set_trace()
    
    assert ndj['reference']['nuc']==ref, "JSON reference doesn't match provided FASTA"

    # Format of vcf seqs: { 'seq1':{4:'A', 7:'-'}, 'seq2':{100:'C'} }, (0-based)
    v_var = vcf['sequences']
    # The parsed VCF only details changes from the ref, but we need
    # to include the ref alleles to compare to the JSON
    for name, data in v_var.items():
        for pos in vcf['positions']: # 0-based
            if pos not in data:
                data[pos] = ref[pos]

    for mut in ndj['nodes'][T.root.name].get('muts', []):
        (frm,pos,to) = (mut[0], int(mut[1:-1])-1, mut[-1])
        assert ref[pos]==frm, f"JSON: Mutation at root node {mut!r} disagrees with reference {ref[pos]!r}"

    j_var = defaultdict(dict)
    j_positions = set() # 0-based
    for clade in T.find_clades():
        # print("Clade", clade)
        muts = ndj['nodes'][clade.name].get('muts', [])
        for descendant in clade.find_clades(): # includes clade itself
            # print(f"\t{descendant}")
            for mut in muts:
                (frm,pos,to) = (mut[0], int(mut[1:-1])-1, mut[-1])
                j_positions.add(pos)
                j_var[descendant.name][pos]=to
    
    # For variable sites, fill in the nodes above the mutations
    # with the reference base
    for clade in T.find_clades():
        for pos in j_positions:
            if pos not in j_var[clade.name]:
                j_var[clade.name][pos] = ref[pos]

    problems_found = False

    for name in set(list(j_var.keys()) + list(v_var.keys())):
        vp = set(v_var[name].keys())
        jp = set(j_var[name].keys())
        if len(vp - jp):
            print(f"{name} had extra positions in VCF")
            problems_found=True
        if len(jp - vp):
            print(f"{name} had extra positions in JSON")
            problems_found=True
        for pos in jp.union(vp):
            try:
                jallele = j_var[name][pos]
            except KeyError:
                print(f"{name}. VCF had pos {pos+1} but JSON didn't. Allele: {v_var[name][pos]}. JSON mask: {bool(mask[pos])}")
                problems_found=True
                continue
            try:
                vallele = v_var[name][pos]
            except KeyError:
                print(f"{name}. JSON had pos {pos+1} but VCF didn't. Allele: {jallele}")
                problems_found=True
                continue
            if jallele!=vallele:
                print(f"{name} pos {pos+1} VCF: {vallele} JSON: {jallele}")
                problems_found=True
    sys.exit(2 if problems_found else 0)