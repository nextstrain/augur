import matplotlib as mpl
from matplotlib import pyplot as plt
from datetime import datetime

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
                lat_long_db[line['location']] = {
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
