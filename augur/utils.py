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

def read_nodedata(fname, traits=None):
    import json
    if os.path.isfile(fname):
        with open(fname) as jfile:
            nodedata = json.load(jfile)
    else:
        print("ERROR: node data can't be read, file %s not found"%fname)
    if traits and os.path.isfile(traits):
        with open(traits) as jfile:
            trait_data = json.load(jfile)
        for k,v in trait_data.items():
            if k in nodedata["nodes"]:
                nodedata["nodes"][k].update(v)
    return nodedata


def write_json(data, file_name, indent=1):
    import json
    import os

    #in case auspice folder does not exist yet
    if not os.path.exists(os.path.dirname(file_name)):
        try:
            os.makedirs(os.path.dirname(file_name))
        except OSError: #Guard against race condition
            if not os.path.isdir(os.path.dirname(file_name)):
                raise
    try:
        handle = open(file_name, 'w')
    except IOError:
        raise
    else:
        json.dump(data, handle, indent=indent)
        handle.close()
