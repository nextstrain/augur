from Bio import SeqIO
import pandas as pd
from pandas import DataFrame
import os.path
from os.path import isfile
import collections
import numpy as np
import geopy
from geopy.geocoders import Nominatim
from datetime import datetime
from matplotlib import cm
from .utils import ambiguous_date_to_date_range, read_lat_longs

forbidden_characters_noFasta = [(' ',''), ('(','_'),(')','_'),(':','_'),(',','_'),(';','_'),('\\','_')]
forbidden_characters = [(' ',''), ('(','_'),(')','_'),(':','_'),(',','_'),(';','_'),('\\','_'),('-','_')]
indicators_for_unknown = ['', '?', 'unknown', 'na', 'nan', '-n/a-', '-']
default_unknown = ''

def determine_standard_columns(columns):
    #force user to keep required columns
    print('Please identify columns containing required meta data. \nThese are the available columns:')
    for i, j in enumerate(columns):
        print(i, j)
    print(i+1, 'none of them')
    take = {}
    for field, explanation in [("strain",  "(matches sequence names)"), ("accession",""),
                               ("date",""), ("country","")]:
        ti = int(input('Which one of your columns is "%s"%s?'%(field, explanation)))
        print()
        if ti<len(columns):
            take[field] = columns[ti]
    return take


def determine_columns_to_keep(remaining_columns):
    print('\n\nRegion information will be generated within this script.')
    print('If you have region information in the data it will be in conflict with regions generated here.')
    print('One solution is to not take column "region".')
    print('These are the available columns not already taken:')
    for j in remaining_columns:
        print('- ',j)

    columns_to_keep = {}
    for j in remaining_columns:
        print("Take column \'{}\'?".format(j))
        c = input("[y] or n: ")

        if c.lower() != 'n':
            columns_to_keep[j] = j

    return columns_to_keep


#Can these be replaced with functions from util?
def numerical_date(date):
    from datetime import datetime
    days_in_year = date.toordinal()- datetime(year=date.year, month=1, day=1).date().toordinal()
    return date.year + days_in_year/365.25

#Can this be replaced with functions from util?
def parse_date(datein, fmt, minmax):
    from datetime import datetime
    import numpy as np
    try:
        if 'XX' in datein:
            min_date, max_date = ambiguous_date_to_date_range(datein, fmt, minmax)
            n_date = np.array((numerical_date(min_date), numerical_date(max_date)))
        else:
            tmp = datetime.strptime(datein, fmt).date()
            n_date = numerical_date(tmp)
    except:
        print("Can't parse ",datein)
        n_date=None

    return n_date

# ...
# changing date format to %Y-%m-%d
# ...
def adjust_date(strains):
    import random
    dropped_strains = []
    print("sample of date specs:")
    for x in random.sample(list(strains.values()), 10):
        print(x['date'])

    print('Example date format: %m-%d-%Y')
    date_format = input("Enter date format (partial dates will be parsed separately): ")
    separator = date_format.split('%')[1][1:2]
    minmax = [1700, 2100]  #set this to be something better

    for s, v in strains.items():
        #
        if len(v['date'])==4: # check whether this works as a number like 2005
            v['date'] = v['date'] + '-XX-XX'
            if parse_date(v['date'], '%Y-%m-%d', minmax)[0] == None:
                dropped_strains.append(s)

        elif len(v['date']) == 7: # check whether this can be parsed as 2005-04 or 04-2004
            if 'm' in date_format.split('Y')[0]:
                v['date'] = v['date'][3:] + '-' + v['date'][:2] + '-XX'
            else:
                v['date'] = v['date'][:4] + '-' + v['date'][5:7] + '-XX'

            if parse_date(v['date'], '%Y-%m-%d', minmax)[0] == None:
                dropped_strains.append(s)

        elif 'X' in v['date']:
            dec_date = parse_date(v['date'], date_format, minmax)
            if dec_date[0] == None:
                dropped_strains.append(s)
            else:
                part1 = v['date'].split(separator)[0]
                part2 = v['date'].split(separator)[1]
                part3 = v['date'].split(separator)[2]
                if len(part1) == 4:
                    year = part1
                    if part2 != 'XX':
                        month = part2
                    else:
                        month = part3
                else:
                    if len(part3) == 4:
                        year = part3
                        if part1 != 'XX':
                            month = part1
                        else:
                            month = part2

                v['date'] = year + '-' + month + '-XX'

        else: #looks like a proper date, parse and recast into ISO format
            dec_date = parse_date(v['date'], date_format, minmax)
            if dec_date == None:
                dropped_strains.append(s)
            else:
                v['date'] = datetime.strptime(v['date'], date_format).strftime("%Y-%m-%d")

    return(strains, dropped_strains)


def read_region_info(region_file=None):
    region_map = {}
    if region_file:
        with open(args.geo_regions, mode='r') as f:
            for line in f:
                entries = line.strip().split('\t')
                region_map[entries[0]] = entries[1]
    else:
        from pkg_resources import resource_stream
        with resource_stream(__package__, "data/geo_regions.tsv") as stream:
            for line in stream:
                entries = line.strip().split()
                if entries:
                    region_map[entries[0]] = entries[1]
    return region_map


def get_synonyms(syn_path):
    synonyms = {}
    #read in in more human-readable format
    with open(syn_path, mode='rt') as f:
        for line in f:
            vals = line.strip().split('\t')
            key = vals[0]
            syms = vals[1:]
            synonyms[key] = syms

    #reorganize to be easier to use in code
    syns = {}
    for key, val in synonyms.items():
        for j in val:
            syns[j] = key

    return syns

# ...
# asks user to match country to region
# ...
def region_questionnaire(sample, region_dict, regions):

    if sample['country'] in region_dict:
        sample['region'] = region_dict[country]
    else:
        print("Available regions are:")
        for i, j in enumerate(regions):
            print(i, j)

        print("Unknown location:", sample['country'])
        reg_no = input("Type the region number to specify which region the country belongs to: ")
        sample['region'] = regions[int(reg_no)]

        region_dict[sample['country']] = sample['region']


def geo_color_generator(countries_current, countries_by_region, region_dict, region_colors):

    colors = []

    for place in countries_current:
        html_col = []
        region = region_dict[place]
        area_diversity = 4 * len(countries_by_region[region])
        dev = np.random.normal(0, area_diversity, 3)

        html_col_r = int(region_colors[region][0][1:3], 16) + int(dev[0])
        html_col_g = int(region_colors[region][0][3:5], 16) + int(dev[1])
        html_col_b = int(region_colors[region][0][5:7], 16) + int(dev[2])
        if html_col_r < 0:
            html_col_r = 00
        if html_col_r > 255:
            html_col_r = 255
        if html_col_g < 0:
            html_col_g = 00
        if html_col_g > 255:
            html_col_g = 255
        if html_col_b < 0:
            html_col_b = 00
        if html_col_b > 255:
            html_col_b = 255
        html_col.append(html_col_r)
        html_col.append(html_col_g)
        html_col.append(html_col_b)
        c_c = "#"+"".join("%02hx"%x for x in html_col)
        colors.append((place, c_c))

    return colors


def write_synonyms(synonyms, syn_path):
    #convert back to right format
    synList = {}
    for key, val in synonyms.items():
        if val in synList.keys():
            synList[val].append(key)
        else:
            synList[val] = [key]

    #put together to write out
    syn_write = []
    for key, val in synList.items():
        line = [key] + val
        fin_lin = "\t".join(line)
        syn_write.append(fin_lin)

    with open(syn_path, 'w') as the_file:
        the_file.write("\n".join(syn_write))


def locate_place(sample, coor, geolocator, region_dict, regions, synonyms):
    place_name = sample['country']
    location = geolocator.geocode((place_name).replace('_',' '), language='en')
    if str(location) == 'None':
        print("The location as formatted couldn't be found in GeoPy")
        answer = 'n'
    else:
        new_loc = (location.address).lower().replace(' ','_')
        print("the unknown location is: ", place_name)
        print("GeoPy suggests using:", new_loc)
        answer = input("Is this the right country name? y or n: ")

    if answer.lower() == 'y':
        sample['country'] = new_loc
        if new_loc not in coor:
            lon = location.longitude
            lat = location.latitude
            coor[new_loc] = [new_loc, 'XX', lat, lon]

        region_questionnaire(sample, region_dict, regions)

    else:
        # Here we let the user try to type something else to find in synonyms or geoPy
        redo = True
        while redo == True:
            print(place_name)
            new_loc = (input("Type the correct country name or 'n' if unknown: ")).lower().replace(' ','_')
            if new_loc in ['n', 'na']:
                indicators_for_unknown.append(place_name)
                sample['country'] = default_unknown
                sample['region'] = default_unknown

            else:
                new_loc = new_loc.lower().replace(' ','_')
                if new_loc in synonyms.keys():
                    new_loc = synonyms[new_loc]

                sample['country'] = new_loc

                if new_loc in region_dict: #what was typed, or synonym of what was typed, is in geo
                    sample['region'] = region_dict[new_loc]
                    redo = False
                else:
                    if new_loc != 'n' and new_loc in coor:
                        region_questionnaire(sample, region_dict, regions)
                        redo = False

                    elif new_loc != 'n' and new_loc not in coor:
                        location = geolocator.geocode((new_loc).replace('_',' '), language='en')
                        print('would you like to: \n1. Use the coordinates for this location?', location.address)
                        print('2. Enter coordinates manually in the file after this run')
                        print('3. Type the country again')
                        loc_ok = input('1, 2, or 3: ')
                        if loc_ok == '1':
                            lon = location.longitude
                            lat = location.latitude
                            coor[new_loc] = [new_loc, 'XX', lat, lon]
                            redo = False
                        elif loc_ok == '2':
                            print('update the new longitude and latitude manually in the geo_long_lat file')
                            redo = False
                        #else 3, allow user to type again

                        if loc_ok != '3':
                            print('What name would you like to use?')
                            print('1: ', new_loc)
                            print('2: ', (location.address).lower().replace(' ','_'))
                            name_loc = input('1 or 2: ')
                            if name_loc == '2':
                                sample['country'] = (location.address).lower().replace(' ','_')
                                new_loc = sample['country']
                            region_questionnaire(sample, region_dict, regions)



def run(args):

    geolocator = Nominatim()

    #Read in meta-data provided by user
    orig_meta = pd.read_csv(args.input, sep="\t" if args.input.endswith('.tsv') else r'\s*,\s*')

    #Get user to identify columns, select those to keep
    column_map = determine_standard_columns(orig_meta.columns)

    #Pick other columns to keep
    column_map.update(determine_columns_to_keep([c for c in orig_meta.columns if c not in column_map.values()]))

    ## properties section would go here - leaving out for now.
    ##### (this standardizes your host names, etc)
    #####

    new_meta = orig_meta.loc[:, column_map.values()].rename(column_map)
    new_meta['region'] = 0
    samples = {}
    for ri, r in new_meta.iterrows():
        if r.strain not in samples:
            samples[r.strain] = r.to_dict()
        else:
            print("dropped duplicate strain:", r.strain)

    # Removing forbidden characters from 'strain' goes here
    # But unless you also modify the Fasta/VCF this can cause more trouble than its worth!
    # Add later but only when can modify all files at once....
    ######
    ######

    #Change date format
    '''
    (note from sanda): This is still not so robust, works for stuff like %Y/%m/%d, %m_%d_%Y....
    Doesn't work for stuff like: %Y///%m///%d or 4320319
    For unrecognized formats replace the funct adjust_date by some other one that returns the date as %Y-%m-%d
    '''
    samples, dropped_strains = adjust_date(samples)
    samples = {k:v for k,v in samples.items() if k not in dropped_strains}

    # ...
    # filter meta (drop duplicates and the ones with unparsable dates)
    # ...
    #TODO: what is 'unparsable dates'? Are those with no dates dropped?
    '''
    (note from sanda): this can have additional filtering criteria, just append to drop_tsv
    perhaps import and append indexes from dropped_strains.txt if you need them to be dropped as well
    '''

    ##Get region information - file assigns each country a region
    region_dict = read_region_info()
    if args.region_colors:
        region_colors = pd.read_csv(args.region_colors, sep='\t')
    else:
        region_colors = {}

    ##Get country synonyms. We store in one format (more human friendly), and 'use' in another
    if args.synonyms:
        synonyms = get_synonyms(args.synonyms)
    else:
        synonyms = {}

    # Get lat-longs - needs to be adapted to new format!!
    #coor = pd.read_csv(args.geo_info, sep='\t')
    coordinates = read_lat_longs()  # should use this instead! But needs adapting

    #go through table...
    countries_current = [] #this keeps track for generating colours later...

    for s,v in samples.items():
        wrong_name = ""
        #country = v['country']

        if type(v['country']) == str:
           v['country'] = v['country'].lower().replace(' ','_')
        else:
            v['country'] = default_unknown

        if v['country'] in region_dict:
            v['region'] = region_dict[v['country']]
        elif v['country'] in indicators_for_unknown:
            v['country'] = default_unknown
            v['region'] = default_unknown
        elif v['country'] in synonyms:
            v['country'] = synonyms[v['country']]
            v['region'] = region_dict[v['country']]
        else: # country not the region dict nor in list of synonyms nor unknown
            wrong_name = v['country']
            locate_place(v, coordinates, geolocator,region_dict, region_colors.keys(), synonyms)

            print(wrong_name)
            s_add = input("Was this country misspelled (should it be added to synonyms)? y or n: ")
            if s_add.lower() == 'y':
                synonyms[wrong_name] = v['country']


    countries_current = {v['country'] for v in samples.values()}
    regions_current = {v['region'] for v in samples.values()}
    countries_by_region = defaultdict(list)
    for v in samples.values():
        countries_by_region[v['region']].append(v['country'])

    ## Properties things will go here eventually (host, etc)
    ####
    ####

    # ...
    # generate region/country colors
    # ...
    geo_color = geo_color_generator(countries_current, countries_by_region,
                                    region_dict, region_colors)



    # ...
    # writing out color/meta
    # ...
    mat_reg = region_colors.as_matrix()
    data = np.array([["region"]*len(region_colors.columns)+[""]+["country"]*len(countries_current),
                    list(region_colors.columns)+[""]+countries_current, list(mat_reg[0])+[""]+geo_color])


    #write out colours
    data = data.transpose()
    color_info = pd.DataFrame(data)
    color_info = color_info.to_csv(sep='\t', index=False, header=False) #, header=['property', 'name', 'color'])
    with open(args.output_color, 'w') as outfile:
        outfile.write(color_info)

    #write out meta
    meta_str = new_meta.to_csv(sep='\t', index=False)
    with open(args.output_meta, 'w') as outfile:
       outfile.write(meta_str)

    # ...
    # update files: geo_info/geo_lat_long/properties_info/synonyms
    # ...
    #geo_regions
    (pd.DataFrame.from_dict(data=region_dict, orient='index').to_csv(args.geo_regions, sep='\t', header=False))

    #lat-long
    #TODO: Need to decide if we want to update the augur lat-long, or just a local version?
    #Does this not do anything??
    with open(args.geo_info, 'w') as outfile:
        outfile.write(coor.to_csv(sep='\t', index=False, header=coor.columns))

    #if include hosts/properties someday, do here

    #synonyms
    write_synonyms(synonyms, args.synonyms)

    #filtering and writing fasta would go here if included.
