from Bio import Entrez, SeqIO
from StringIO import StringIO

Entrez.email = "richard.neher@unibas.ch"     # Always tell NCBI who you are

outgroups = {#'H3N2':'U26830',
             #'H3N2_77':'CY113261.1',
             #'H1N1pdm':'AF455680',
             #'Vic':'CY018813',
             #'Yam':'CY019707',
             #'zika':'KX369547.1'
             'yam_pb1':'CY019713',
             'yam_pb2':'CY019714',
             'yam_pa': 'CY019712',
             'yam_ha': 'CY019707',
             'yam_np': 'CY019710',
             'yam_na': 'CY019709',
             'yam_ma': 'CY019708',
             'yam_ns': 'CY019711',

             'vic_pb1':'CY018819',
             'vic_pb2':'CY018820',
             'vic_pa': 'CY018818',
             'vic_ha': 'CY018813',
             'vic_np': 'CY018816',
             'vic_na': 'CY018815',
             'vic_ma': 'CY018814',
             'vic_ns': 'CY018817',

             'h1n1pdm_pb1':'AF455728',
             'h1n1pdm_pb2':'AF455736',
             'h1n1pdm_pa': 'AF455720',
             'h1n1pdm_ha': 'AF455680',
             'h1n1pdm_np': 'AF455704',
             'h1n1pdm_na': 'AF455696',
             'h1n1pdm_ma': 'AF455688',
             'h1n1pdm_ns': 'AF455712',

             'h3n2_pb1':'CY113683',
             'h3n2_pb2':'CY113684',
             'h3n2_pa': 'CY113682',
             'h3n2_ha': 'CY113677',
             'h3n2_np': 'CY113680',
             'h3n2_na': 'CY113679',
             'h3n2_ma': 'CY113678',
             'h3n2_ns': 'CY113681',
             }

if __name__ == "__main__":
    for virus, genbank_id in outgroups.iteritems():
        handle = Entrez.efetch(db="nucleotide", id=genbank_id, rettype="gb")
        seq= SeqIO.read(StringIO(handle.read()), format = 'genbank')
        SeqIO.write(seq, 'flu/metadata/'+virus+'_outgroup.gb', format='genbank')
