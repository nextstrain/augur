import pandas as pd
import numpy as np

CDC_h1n1pdm_antigenic_table = 'metadata/igarashi_H1N1_antigenic_sites.xlsx'

df = pd.read_excel(CDC_h1n1pdm_antigenic_table)

epitope_col = u'Antigenic Site (Caton)'

# there are two extra rows in the file, but this doesn't matter since
# we restrict to length of the total peptide in augur
epi_mask = "".join(['1' if type(x)==unicode else "0" for x in df.loc[:,epitope_col]][1:])
ha1_mask = "".join(['1' if str(x).startswith('HA1') else '0' for x in df.index][1:])
ha1_head_mask = "".join(['1' if str(x).startswith('HA1_globular_head') else '0' for x in df.index][1:])

with open('metadata/ha_masks.tsv', 'a') as ofile:
	ofile.write('canton\t'+epi_mask+'\n')
	ofile.write('ha1_h1n1pdm\t'+ha1_mask+'\n')
	ofile.write('ha1_globular_head_h1n1pdm\t'+ha1_head_mask+'\n')
