from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation
from glob import glob

record_paths = glob('/Users/Sidney/nextstrain/augur/dengue/metadata/*.gb')
proteins = ['C', 'M', 'E', 'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', '2K', 'NS4B', 'NS5']
# refseq_names = {'NC_001477': 'DENV1/NAURUISLAND/REFERENCE/1997', #DENV1
#                 'NC_001474': 'DENV2/THAILAND/REFERENCE/1964', #DENV2
#                 'NC_001475': 'DENV3/SRI_LANKA/REFERENCE/2000', #DENV3
#                 'NC_002640': 'DENV4/NA/REFERENCE/2003' #DENV4
#                }

refseq_names = {'NC_001477': 'DENV1', #DENV1
                'NC_001474': 'DENV2', #DENV2
                'NC_001475': 'DENV3', #DENV3
                'NC_002640': 'DENV4' #DENV4
               }
print 'WARNING: locus field must be edited by hand to match the strain name in vdb. (Too long for BioPython to write).'

def convert_annotations(record_path):
    record = SeqIO.read(record_path, 'genbank')
    record.name = refseq_names[record.id.split('.')[0]]

    seen_proteins = [] # Keep the first instance of a gene name to maintain the longest entry (e.g., M instead of prM)
    for f in record.features:
        if 'gene' in f.qualifiers and 'product' in f.qualifiers:
            protein_product = f.qualifiers['product'][0].split()[-1]
            if protein_product not in seen_proteins:
                f.qualifiers['gene'] = [protein_product]
                seen_proteins.append(protein_product)
    SeqIO.write(record, record_path, 'genbank')

for r in record_paths:
    convert_annotations(r)
