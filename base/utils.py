import matplotlib as mpl
from matplotlib import pyplot as plt

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
