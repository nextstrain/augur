import os
import pandas as pd
from treetime.utils import numeric_date

def ambiguous_date_to_date_range(mydate, fmt):
    from datetime import datetime
    sep = fmt.split('%')[1][-1]
    min_date, max_date = {}, {}
    today = datetime.today().date()

    for val, field  in zip(mydate.split(sep), fmt.split(sep+'%')):
        f = 'year' if 'y' in field.lower() else ('day' if 'd' in field.lower() else 'month')
        if 'XX' in val:
            if f=='year':
                return None, None
            elif f=='month':
                min_date[f]=1
                max_date[f]=12
            elif f=='day':
                min_date[f]=1
                max_date[f]=31
        else:
            min_date[f]=int(val)
            max_date[f]=int(val)
    max_date['day'] = min(max_date['day'], 31 if max_date['month'] in [1,3,5,7,8,10,12]
                                           else 28 if max_date['month']==2 else 30)
    lower_bound = datetime(year=min_date['year'], month=min_date['month'], day=min_date['day']).date()
    upper_bound = datetime(year=max_date['year'], month=max_date['month'], day=max_date['day']).date()
    return (lower_bound, upper_bound if upper_bound<today else today)

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
                try:
                    numerical_dates[k] = numeric_date(datetime.strptime(v, fmt))
                except:
                    numerical_dates[k] = None
    else:
        numerical_dates = {k:float(v) for k,v in dates.items()}

    return numerical_dates

def read_node_data(fname, traits=None, aa_muts=None):
    import json
    if os.path.isfile(fname):
        with open(fname) as jfile:
            node_data = json.load(jfile)

        for more_data in [traits, aa_muts]:
            if more_data and os.path.isfile(more_data):
                with open(more_data) as jfile:
                    tmp_data = json.load(jfile)
                for k,v in tmp_data.items():
                    if k in node_data["nodes"]:
                        node_data["nodes"][k].update(v)
    else:
        print("ERROR: node data can't be read, file %s not found"%fname)
        node_data=None

    return node_data


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


def load_features(reference, feature_names=None):
    #read in appropriately whether GFF or Genbank
    #checks explicitly for GFF otherwise assumes Genbank
    if not os.path.isfile(reference):
        print("ERROR: reference sequence not found. looking for", reference)
        return None

    features = {}
    if '.gff' in reference.lower():
        #looks for 'gene' and 'gene' as best for TB
        from BCBio import GFF
        limit_info = dict( gff_type = ['gene'] )

        with open(reference) as in_handle:
            for rec in GFF.parse(in_handle, limit_info=limit_info):
                for feat in rec.features:
                    if "gene" in feat.qualifiers:
                        fname = feat.qualifiers["gene"][0]
                    else:
                        fname = feat.qualifiers["locus_tag"][0]
                    if feature_names is None or fname in feature_names:
                        features[fname] = feat

            if feature_names is not None:
                for fe in feature_names:
                    if fe not in features:
                        print("Couldn't find gene {} in GFF or GenBank file".format(fe))

    else:
        from Bio import SeqIO
        for feat in SeqIO.read(reference, 'genbank').features:
            if feat.type=='CDS':
                if "locus_tag" in feat.qualifiers:
                    fname = feat.qualifiers["locus_tag"][0]
                    if feature_names is None or fname in feature_names:
                        features[fname] = feat
                elif "gene" in feat.qualifiers:
                    fname = feat.qualifiers["gene"][0]
                    if feature_names is None or fname in feature_names:
                        features[fname] = feat

    return features
