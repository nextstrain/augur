import pandas as pd
import numpy as np

# run this script and append output to metadata/ha_masks.tsv

CDC_h1n1pdm_antigenic_table = 'igarashi_H1N1_antigenic_sites.xlsx'

df = pd.read_excel(CDC_h1n1pdm_antigenic_table)

epitope_col = u'Antigenic Site (Caton)'

# Find all rows for which the amino acid position is defined.
aa_position_defined = ~pd.isnull(df["Amino Acid Position"])

# Convert missing (null) values in the antigenic site column to 0 and present values to 1.
epi_mask = "".join((~pd.isnull(df.loc[aa_position_defined, epitope_col])).astype(int).astype(str))

# Find all sites in HA1, generally, and the globular head, specifically.
ha1_mask = "".join(df[aa_position_defined].index.str.startswith("HA1").astype(int).astype(str))
ha1_head_mask = "".join(df[aa_position_defined].index.str.startswith("HA1_globular_head").astype(int).astype(str))

print('canton\t'+epi_mask)
print('ha1_h1n1pdm\t'+ha1_mask)
print('ha1_globular_head_h1n1pdm\t'+ha1_head_mask)
