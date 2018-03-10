import pandas as pd
import numpy as np

CDC_h1n1pdm_antigenic_table = 'metadata/igarashi_H1N1_antigenic_sites.xlsx'

df = pd.read_excel(CDC_h1n1pdm_antigenic_table)

epitope_col = u'Antigenic Site (Caton)'

# there are two extra rows in the file, hence -2
mask = "".join([ '1' if type(x)==unicode else "0" for x in df.loc[:-2,epitope_col]])

with open('metadata/ha_masks.tsv', 'a') as ofile:
	ofile.write('canton\t'+mask+'\n')
