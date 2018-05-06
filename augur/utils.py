import pandas as pd
from treetime.utils import numeric_date
from base.utils import ambiguous_date_to_date_range

def parse_metadata(fname):
    metadata = pd.read_csv(fname, sep='\t' if fname[-3:]=='tsv' else ',',
                         skipinitialspace=True)
    return metadata

def meta_to_date_dict(metadata, name_col = None, date_col='date', fmt=None):
    possible_name_cols = ['strain', 'name'] if name_col is None else [name_col]
    dates = None
    for nc in possible_name_cols:
        try:
            dates = {row.loc[nc]:row.loc[date_col] for ri,row in metadata.iterrows()}
            break
        except:
            pass

    if dates is None:
        print("ERROR: Can not find meta data column with strain names")
        return None

    if fmt:
        from datetime import datetime
        numerical_dates = {}
        for k,v in dates.items():
            if 'XX' in v:
                numerical_dates[k] = [numeric_date(d) for d in ambiguous_date_to_date_range(v, fmt)]
            else:
                numerical_dates[k] = numeric_date(datetime.strptime(v, fmt))
    else:
        numerical_dates = {k:float(v) for k,v in dates.items()}

    return numerical_dates
