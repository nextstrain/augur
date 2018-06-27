import os, json
import pandas as pd
from treetime.utils import numeric_date
from collections import defaultdict

def myopen(fname, mode):
    if fname.endswith('.gz'):
        import gzip
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)

def ambiguous_date_to_date_range(mydate, fmt, min_max_year=None):
    from datetime import datetime
    sep = fmt.split('%')[1][-1]
    min_date, max_date = {}, {}
    today = datetime.today().date()

    for val, field  in zip(mydate.split(sep), fmt.split(sep+'%')):
        f = 'year' if 'y' in field.lower() else ('day' if 'd' in field.lower() else 'month')
        if 'XX' in val:
            if f=='year':
                if min_max_year:
                    min_date[f]=min_max_year[0]
                    if len(min_max_year)>1:
                        max_date[f]=min_max_year[1]
                    elif len(min_max_year)==1:
                        max_date[f]=4000 #will be replaced by 'today' below.
                else:
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
                             skipinitialspace=True).fillna('')
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


def get_numerical_dates(meta_dict, name_col = None, date_col='date', fmt=None, min_max_year=None):
    if fmt:
        from datetime import datetime
        numerical_dates = {}
        for k,m in meta_dict.items():
            v = m[date_col]
            if type(v)!=str:
                print("WARNING: %s has an invalid data string:"%k,v)
                continue
            elif 'XX' in v:
                ambig_date = ambiguous_date_to_date_range(v, fmt, min_max_year)
                if ambig_date is None or None in ambig_date:
                    numerical_dates[k] = [None, None] #don't send to numeric_date or will be set to today
                else:
                    numerical_dates[k] = [numeric_date(d) for d in ambig_date]
            else:
                try:
                    numerical_dates[k] = numeric_date(datetime.strptime(v, fmt))
                except:
                    numerical_dates[k] = None
    else:
        numerical_dates = {k:float(v) for k,v in meta_dict.items()}

    return numerical_dates

def read_node_data(fnames):
    """parse the "nodes" field of the given JSONs and join the data together"""
    if type(fnames) is str:
        fnames = [fnames]
    node_data = {"nodes": {}}
    for fname in fnames:
        if os.path.isfile(fname):
            with open(fname) as jfile:
                tmp_data = json.load(jfile)
            for k,v in tmp_data.items():
                if k=="nodes":
                    for n,nv in v.items():
                        if n in node_data["nodes"]:
                            node_data["nodes"][n].update(nv)
                        else:
                            node_data["nodes"][n] = nv
                elif k in node_data:
                    # will this recurse into nested dicts?!?!
                    node_data[k].update(v)
                else:
                    node_data[k]=v
        else:
            print("ERROR: node_data JSON file %s not found. Attempting to proceed without it."%fname)
    return node_data


# put the data saved in node_data json back onto the tree
def attach_tree_meta_data(T, node_meta):
    for n in T.find_clades(order='preorder'):
        if n.name not in node_meta:
            print("ERROR: keys in tree and node meta data don't match. Node %s is missing"%n.name)
            continue

        n.attr={}
        n.aa_muts={}
        for field, val in node_meta[n.name].items():
            n.__setattr__(field, val)
            if field=='mutations':
                n.muts = val
            if field=='numdate':
                n.num_date = val


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
        try:
            from BCBio import GFF #Package name is confusing - tell user exactly what they need!
        except ImportError:
            print("ERROR: Package BCBio.GFF not found! Please install using \'pip install bcbio-gff\' before re-running.")
            return -1
        limit_info = dict( gff_type = ['gene'] )

        with open(reference) as in_handle:
            for rec in GFF.parse(in_handle, limit_info=limit_info):
                for feat in rec.features:
                    if feature_names is not None: #check both tags; user may have used either
                        if "gene" in feat.qualifiers and feat.qualifiers["gene"][0] in feature_names:
                            fname = feat.qualifiers["gene"][0]
                        elif "locus_tag" in feat.qualifiers and feat.qualifiers["locus_tag"][0] in feature_names:
                            fname = feat.qualifiers["locus_tag"][0]
                        else:
                            fname = None
                    else:
                        if "gene" in feat.qualifiers:
                            fname = feat.qualifiers["gene"][0]
                        else:
                            fname = feat.qualifiers["locus_tag"][0]
                    if fname:
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
            elif feat.type=='source': #read 'nuc' as well for annotations - need start/end of whole!
                features['nuc'] = feat

    return features

def read_config(fname):
    if fname and os.path.isfile(fname):
        with open(fname) as ifile:
            config = json.load(ifile)
    else:
        print("ERROR: config file %s not found."%fname)
        config = defaultdict(dict)

    return config

def read_lat_longs(defaults, overrides=None):
    coordinates = {}
    if defaults:
        if os.path.isfile(defaults):
            with open(defaults) as ifile:
                for line in ifile:
                    if line.startswith('#'): continue
                    fields = line.strip().split()
                    if len(fields) == 4:
                        geo_field, loc = fields[0], fields[1]
                        lat, long = float(fields[2]), float(fields[3])
                        coordinates[(geo_field, loc)] = {
                            "latitude": lat,
                            "longitude": long
                        }
        else:
            print("ERROR: default lat/long file %s not found." % defaults)

    if overrides:
        if os.path.isfile(overrides):
            with open(overrides) as ifile:
                for line in ifile:
                    if line.startswith('#'): continue
                    fields = line.strip().split()
                    if len(fields) == 4:
                        geo_field, loc = fields[0], fields[1]
                        lat, long = float(fields[2]), float(fields[3])
                        coordinates[(geo_field, loc)] = {
                            "latitude": lat,
                            "longitude": long
                        }
        else:
            print("WARNING: input lat/long file %s not found." % overrides)
    return coordinates

def read_colors(defaults, overrides=None):
    colors = {}
    if defaults:
        if os.path.isfile(defaults):
            with open(defaults) as fh:
                for line in fh:
                    if line.startswith('#'): continue
                    fields = line.strip().split()
                    if len(fields) == 3:
                        trait, trait_value, hex_code = fields[0], fields[1], fields[2]
                        colors[(trait, trait_value)] = hex_code
        else:
            print("WARNING: Couldn't open color definitions file {}.".format(defaults))
    if overrides:
        if os.path.isfile(overrides):
            with open(overrides) as fh:
                for line in fh:
                    if line.startswith('#'): continue
                    fields = line.strip().split()
                    if len(fields) == 3:
                        trait, trait_value, hex_code = fields[0], fields[1], fields[2]
                        colors[(trait, trait_value)] = hex_code
        else:
            print("WARNING: Couldn't open color definitions file {}.".format(overrides))
    color_map = defaultdict(list)
    for (trait, trait_value), hex_code in colors.items():
        color_map[trait].append((trait_value, hex_code))
    return color_map

def write_VCF_translation(prot_dict, vcf_file_name, ref_file_name):
    """
    Writes out a VCF-style file (which seems to be minimally handleable
    by vcftools and pyvcf) of the AA differences between sequences and the reference.
    This is a similar format created/used by read_in_vcf except that there is one
    of these dicts (with sequences, reference, positions) for EACH gene.

    Also writes out a fasta of the reference alignment.

    EBH 12 Dec 2017
    """
    import numpy as np

    #for the header
    seqNames = list(prot_dict[list(prot_dict.keys())[0]]['sequences'].keys())

    #prepare the header of the VCF & write out
    header=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+seqNames
    with open(vcf_file_name, 'w') as the_file:
        the_file.write( "##fileformat=VCFv4.2\n"+
                        "##source=NextStrain_Protein_Translation\n"+
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        the_file.write("\t".join(header)+"\n")

    refWrite = []
    vcfWrite = []

    #go through for every gene/protein
    for fname, prot in prot_dict.items():
        sequences = prot['sequences']
        ref = prot['reference']
        positions = prot['positions']

        #write out the reference fasta
        refWrite.append(">"+fname)
        refWrite.append(ref)

        #go through every variable position
        #There are no deletions here, so it's simpler than for VCF nuc sequenes!
        for pi in positions:
            pos = pi+1 #change numbering to match VCF not python
            refb = ref[pi] #reference base at this position

            #try/except is (much) faster than list comprehension!
            pattern = []
            for k,v in sequences.items():
                try:
                    pattern.append(sequences[k][pi])
                except KeyError:
                    pattern.append('.')
            pattern = np.array(pattern)

            #get the list of ALTs - minus any '.'!
            uniques = np.unique(pattern)
            uniques = uniques[np.where(uniques!='.')]

            #Convert bases to the number that matches the ALT
            j=1
            for u in uniques:
                pattern[np.where(pattern==u)[0]] = str(j)
                j+=1
            #Now convert these calls to #/# (VCF format)
            calls = [ j+"/"+j if j!='.' else '.' for j in pattern ]
            if len(uniques)==0:
                print ("UNEXPECTED ERROR WHILE CONVERTING TO VCF AT POSITION {}".format(str(pi)))
                break

            #put it all together and write it out
            output = [fname, str(pos), ".", refb, ",".join(uniques), ".", "PASS", ".", "GT"] + calls

            vcfWrite.append("\t".join(output))

    #write it all out
    with open(ref_file_name, 'w') as the_file:
        the_file.write("\n".join(refWrite))

    with open(vcf_file_name, 'a') as the_file:
        the_file.write("\n".join(vcfWrite))

    if vcf_file_name.lower().endswith('.gz'):
        import os
        #must temporarily remove .gz ending, or gzip won't zip it!
        os.rename(vcf_file_name, vcf_file_name[:-3])
        call = ["gzip", vcf_file_name[:-3]]
        os.system(" ".join(call))
