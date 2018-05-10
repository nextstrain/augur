import os
import pandas as pd
from treetime.utils import numeric_date
from base.utils import ambiguous_date_to_date_range

def read_metadata(fname):
    if os.path.isfile(fname):
        metadata = pd.read_csv(fname, sep='\t' if fname[-3:]=='tsv' else ',',
                             skipinitialspace=True)
        meta_dict = {}
        for ii, val in metadata.iterrows():
            if hasattr(val, "strain"):
                meta_dict[val.strain] = val.to_dict()
            elif hasattr(val, "name"):
                meta_dict[val.name] = val.to_dict()
            else:
                print("ERROR: meta data file needs 'name' or 'strain' column")

        return meta_dict, list(metadata.columns)
    else:
        print("ERROR: file with states does not exist")
        return {}, []


def get_numerical_dates(meta_dict, name_col = None, date_col='date', fmt=None):
    if fmt:
        from datetime import datetime
        numerical_dates = {}
        for k,m in meta_dict.items():
            v = m[date_col]
            if 'XX' in v:
                numerical_dates[k] = [numeric_date(d) for d in ambiguous_date_to_date_range(v, fmt)]
            else:
                numerical_dates[k] = numeric_date(datetime.strptime(v, fmt))
    else:
        numerical_dates = {k:float(v) for k,v in dates.items()}

    return numerical_dates
