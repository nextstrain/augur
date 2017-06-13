import matplotlib as mpl
from matplotlib import pyplot as plt
from datetime import datetime
from base.io_util import tree_to_json
import numpy as np

def generate_cmap(data, discrete):
    '''
    data is set or list (e.g. of countries)
    discrete is bool
    returns dict of (e.g.) country -> hex
    '''
    norm = mpl.colors.Normalize(0, len(data) - 1)
    cmap = mpl.cm.get_cmap("viridis")
    if discrete:
        if len(data) <= 10:
            cmap = mpl.cm.get_cmap("Vega10")
        elif len(data) <= 20:
            cmap = mpl.cm.get_cmap("Vega20")

    ret = {}
    for idx, val in enumerate(list(data)):
        ret[val] = mpl.colors.to_hex(cmap(norm(idx)))
    return ret

def define_latitude_longitude(lat_long_defs, log):
    import csv
    # lat_long_defs is from the config file
    # get the latitude and longitudes that were already determined
    lat_long_db = {}
    if isinstance(lat_long_defs, (list, tuple)):
        paths = lat_long_defs
    else:
        paths = [lat_long_defs]
    for path in paths:
        try:
            file = open(path, 'r')
        except IOError:
            log.warn("Couldn't open lat/long file {}".format(path))
            continue
        reader = csv.DictReader(filter(lambda row: row[0]!='#', file), delimiter='\t')		# list of dicts
        for line in reader:
            try:
                lat_long_db[line['location'].lower()] = {
                    'latitude': float(line['latitude']),
                    'longitude': float(line['longitude'])
                    }
            except:
                log.warn("Error reading line in lat/long file {}. Line: {}".format(path, line))
                continue
        file.close()
    return lat_long_db


def fix_names(n):
    return n.replace(" ","_").replace("(",'_').replace(")",'_').replace("'",'_').replace(":",'_')
#
# def calc_af(aln, alpha):
#     aln_array = np.array(aln)
#     af = np.zeros((len(alpha), aln_array.shape[1]))
#     for ai, state in enumerate(alpha):
#         af[ai] += (aln_array==state).mean(axis=0)
#     af[-1] = 1.0 - af[:-1].sum(axis=0)
#     return af

def num_date(date):
    days_in_year = date.toordinal()- datetime(year=date.year, month=1, day=1).date().toordinal()
    return date.year + days_in_year/365.25

def ambiguous_date_to_date_range(mydate, fmt):
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

def save_as_nexus(tree, fname):
    def format_string_attr(node, key):
        return "{}=\"{}\"".format(key, node["attr"][key])

    def stringify_node(node, prev_div, terminal):
        if terminal:
            taxa.append(node["strain"])
        extra = [format_string_attr(node, x) for x in attrs_to_write]
        return "{}[&{}]:{}".format(len(taxa) if terminal else "", ",".join(extra), float(node["attr"]["div"]) - float(prev_div))

    def tree_walk(node, prev_div):
        if "children" in node:
            subtrees = ",".join([tree_walk(child, node["attr"]["div"]) for child in node["children"]])
            return "({}){}".format(subtrees,stringify_node(node, prev_div, False))
        else:
            return stringify_node(node, prev_div, True)

    def nexus(tree):
        nex = []
        nex += ["#NEXUS", ""]
        nex += ["Begin taxa;", "\tDimensions ntax={};".format(len(taxa)), "\tTaxlabels"]
        nex += ["\t\t"+name for name in taxa]
        nex += ["\t\t;", "End;", ""]
        nex += ["Begin trees;", "\tTranslate"]
        nex += ["\t\t{} {},".format(idx+1, name) for idx, name in enumerate(taxa)]
        nex[-1] = nex[-1][:-1] # remove the final comma
        nex += ["\t\t;", "tree TREE1 = [&R] "+tree+";"]
        nex += ["End;"]
        return nex

    print("Saving to nexus tree {}".format(fname))
    attrs_to_write = tree.root.attr.keys()
    taxa = [] # the id of a name is its (0 based) idx + 1, populated by tree_walk()
    json = tree_to_json(tree.root , ['clade', 'attr']);
    nex = nexus(tree_walk(json, 0))
    with open(fname, 'w') as f:
        f.write("\n".join(nex))


"""
given a date string and a (list of) possible formats, return a list of length three consisting of:
(1) the original string passed in (normally saved as raw_date)
(2) the numerical date - either a float or an array of 2 floats (normally saved as num_date)
(3) the date in a datetime object (normally saved as date)
"""
def parse_date(datein, fmts):
    ret = [datein, None, None]
    for fmt in fmts:
        try:
            if 'XX' in datein:
                min_date, max_date = ambiguous_date_to_date_range(datein, fmt)
                ret[1] = np.array((num_date(min_date), num_date(max_date)))
                ret[2] = min_date
            else:
                if callable(fmt):
                    tmp = fmt(datein)
                else:
                    try:
                        tmp = datetime.strptime(datein, fmt).date()
                    except:
                        tmp = datein
                ret[1] = num_date(tmp)
                ret[2] = tmp
                break
        except:
            continue
    return ret
