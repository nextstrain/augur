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
from io import TextIOWrapper
from .utils import ambiguous_date_to_date_range, read_lat_longs
from collections import defaultdict

forbidden_characters_noFasta = [(' ',''), ('(','_'),(')','_'),(':','_'),(',','_'),(';','_'),('\\','_')]
forbidden_characters = [(' ',''), ('(','_'),(')','_'),(':','_'),(',','_'),(';','_'),('\\','_'),('-','_')]
indicators_for_unknown = ['', '?', 'unknown', 'na', 'nan', '-n/a-', '-']
default_unknown = ''

def determine_standard_columns(columns, country):
    #force user to keep required columns
    print('\n\nPlease identify columns containing required meta data. \nThese are the available columns:')
    for i, j in enumerate(columns):
        print(i, j)
    print(i+1, 'none of them')
    take = {}
    last_column = ("country", "") if country else ("location", "")
    standard_columns = [("strain",  " (matches sequence names)"), ("accession",""),
                               ("date",""), last_column]

    for field, explanation in standard_columns:
        ti = int(input('Which one of your columns is "%s"%s? '%(field, explanation)))
        print()
        if ti<len(columns):
            take[columns[ti]] = field
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
# TODO: check if it's %Y-%m-%d already and don't ask questions if it is...
def adjust_date(strains):
    import random
    dropped_strains = []
    print("\nSample of date format:")
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
            with TextIOWrapper(stream, "utf-8") as defaults:
                for line in defaults:
                    entries = line.strip().split()
                    if entries:
                        region_map[entries[0]] = entries[1]
    return region_map

def read_location_info(location_file):
    location_country = {}
    if location_file:
        with open(location_file, mode='r') as f:
            for line in f:
                entries = line.strip().split('\t')
                location_country[entries[0]] = entries[1]
    return location_country

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
        sample['region'] = region_dict[sample['country']]
    else:
        print("\nWhat region does %s belong to? Available regions are:"%sample['country'])
        for i, j in enumerate(regions):
            print(i, j)

        #print("Unknown location:", sample['country'])
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
        #colors.append((place, c_c)) # how currently written out, won't work with place already attached..
        colors.append(c_c)

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

def get_field_to_use(location):
    field_to_use = 'city'
    if 'city' not in location.raw['address']:
        print("'City' not found in the GeoPy results. Available fields:")
        avail_fields = {}
        for i, j in enumerate(location.raw['address']):
            print(i, j, "(",location.raw['address'][j], ")")
            avail_fields[i] = j
        field_no = input("Type what to use for 'city': ")
        field_to_use = avail_fields[int(field_no)]
    return field_to_use

def add_to_coor(location, name, place_type, coor):
    new_key = (place_type, name)
    if new_key not in coor:
        lon = location.longitude
        lat = location.latitude
        coor[new_key] = {'latitude': lat, 'longitude': lon}


def locate_place(sample, coor, geolocator, region_dict, regions, synonyms, country):
    loc_name = 'country' if country else 'location'
    place_name = sample[loc_name]
    location = geolocator.geocode((place_name).replace('_',' '), language='en', addressdetails=True)

    if str(location) == 'None':
        print("The location as formatted couldn't be found in GeoPy")
        answer = 'n'
    else:
        print("\nAn unknown location has been found: ", place_name)
        if country:
            new_loc = (location.address).lower().replace(' ','_')
            print("GeoPy suggests using:", new_loc)
            answer = input("Is this the right country name, and the name you want to use? y or n: ")
        else:
            field_to_use = get_field_to_use(location)
            new_loc = location.raw['address'][field_to_use].lower().replace(' ', '_')
            countryVal = location.raw['address']['country'].lower()
            print("GeoPy suggests using: %s with country: %s"%(new_loc, countryVal))
            answer = input("Is this the right location and country, and the location name you'd like to use? y or n: ")

    if answer.lower() == 'y':
        sample['country'] = new_loc if country else countryVal
        if not country:
            sample['location'] = new_loc
        #add to coordinates if not there
        add_to_coor(location, new_loc, loc_name, coor)
        if not country: #add country if not country
            country_latlong = geolocator.geocode(countryVal, language='en', addressdetails=True)
            add_to_coor(country_latlong, countryVal, 'country', coor)

        region_questionnaire(sample, region_dict, regions)

    else:
        # Here we let the user try to type something else to find in synonyms or geoPy
        redo = True
        while redo == True:
            # print(place_name)
            new_loc = (input("Type the correct name or 'n' if unknown: ")).lower().replace(' ','_')
            if new_loc in ['n', 'na']:
                indicators_for_unknown.append(place_name)
                if not country:
                    sample['location'] = default_unknown
                sample['country'] = default_unknown
                sample['region'] = default_unknown

            else:
                new_loc = new_loc.lower().replace(' ','_')
                if new_loc in synonyms.keys():
                    new_loc = synonyms[new_loc]

                sample[loc_name] = new_loc

                if new_loc in region_dict: #what was typed, or synonym of what was typed, is in geo
                    sample['region'] = region_dict[new_loc]
                    redo = False
                else:
                    if new_loc != 'n' and (loc_name, new_loc) in coor:
                        region_questionnaire(sample, region_dict, regions)
                        redo = False

                    elif new_loc != 'n' and (loc_name, new_loc) not in coor:
                        location = geolocator.geocode((new_loc).replace('_',' '), language='en', addressdetails=True)
                        coord_text = "coordinates" if country else "coordinates/country"
                        print('Would you like to: \n1. Use the %s for this location?'%(coord_text), location.address)
                        print('2. Enter %s manually in the file after this run'%(coord_text))
                        print('3. Type the location again')
                        loc_ok = input('1, 2, or 3: ')
                        if loc_ok == '1':
                            redo = False
                        elif loc_ok == '2':
                            if not country:
                                new_country = input("Please type the country you'd like to use: ")
                                sample['country'] = new_country
                            region_questionnaire(sample, region_dict, regions)
                            print('Please update the new longitude and latitude manually in the geo_long_lat file')
                            redo = False
                        #else 3, allow user to type again

                        if loc_ok != '3':
                            new_name = (location.address).lower().replace(' ','_')
                            if not country:
                                field_to_use = get_field_to_use(location)
                                new_name = location.raw['address'][field_to_use].lower().replace(' ','_')
                            print('What name would you like to use?')
                            print('1: ', new_loc)
                            print('2: ', new_name)
                            name_loc = input('1 or 2: ')
                            if name_loc == '2':
                                sample[loc_name] = new_name
                                new_loc = sample[loc_name]

                            #add to coordinates, using the found location
                            add_to_coor(location, new_loc, loc_name, coor)
                            if not country:
                                countryVal = location.raw['address']['country'].lower()
                                country_latlong = geolocator.geocode(countryVal, language='en', addressdetails=True)
                                add_to_coor(country_latlong, countryVal, 'country', coor)
                                #and add contry info
                                sample['country'] = location.raw['address']['country'].lower()

                            region_questionnaire(sample, region_dict, regions)



def run(args):

    geolocator = Nominatim()

    #Read in meta-data provided by user
    orig_meta = pd.read_csv(args.input, sep="\t" if args.input.endswith('.tsv') else r'\s*,\s*')

    country = True
    #See if they want to do by Country or smaller location
    print('Would you like to identify geography by country or more specific location?')
    print('1. Country')
    print('2. More specific (city, county, state)')
    ti = int(input('Type your choice (default 1): '))
    if ti==2:
        country = False

    #Get user to identify columns, select those to keep
    column_map = determine_standard_columns(orig_meta.columns, country)

    #Pick other columns to keep
    columns_to_keep = determine_columns_to_keep([c for c in orig_meta.columns if c not in column_map.keys()])
    column_map.update(columns_to_keep)

    ## properties section would go here - leaving out for now.
    ##### (this standardizes your host names, etc)
    #####

    new_meta = orig_meta.loc[:, column_map.keys()].rename(index=str, columns=column_map)
    new_meta['region'] = 0
    if not country:
        new_meta['country'] = 0
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

    #if available, read in file linking location to country
    if not country:
        location_country = {}
        if args.location_country:
            location_country = read_location_info(args.location_country)

    # Get lat-longs - needs to be adapted to new format!!
    #coor = pd.read_csv(args.geo_info, sep='\t')
    coordinates = read_lat_longs()  # should use this instead! But needs adapting
    orig_coordinate_length = len(coordinates)

    #go through table...
    countries_current = [] #this keeps track for generating colours later...

    for s,v in samples.items():
        wrong_name = ""
        loc_name = 'country' if country else 'location'
        #country = v['country']

        if type(v[loc_name]) == str:
           v[loc_name] = v[loc_name].lower().replace(' ','_')
        else:
            v[loc_name] = default_unknown

        if v[loc_name] in region_dict:
            v['region'] = region_dict[v[loc_name]]
        elif v[loc_name] in indicators_for_unknown:
            v[loc_name] = default_unknown
            v[loc_name] = default_unknown
        elif not country and v[loc_name] in location_country:
            v['country'] = location_country[v[loc_name]]
            v['region'] = region_dict[v['country']]
        elif v[loc_name] in synonyms:
            v[loc_name] = synonyms[v[loc_name]] #if country else location_country[synonyms[v[loc_name]]]
            if not country:
                v['country'] =  location_country[v[loc_name]]
            v['region'] = region_dict[v[loc_name]] if country else region_dict[v['country']]
        else: # country not the region dict nor in list of synonyms nor unknown
            wrong_name = v[loc_name]
            locate_place(v, coordinates, geolocator,region_dict, region_colors.keys(), synonyms, country)

            #print(wrong_name)
            if wrong_name != v[loc_name]: #don't offer to add to synonyms if identical...
                s_add = input("Was this location (%s) misspelled (should it be added to synonyms)? y or n: " % wrong_name)
                if s_add.lower() == 'y':
                    synonyms[wrong_name] = v[loc_name]
            if not country:
                location_country[v[loc_name]] = v['country']

    countries_current = {v['country'] for v in samples.values() if v['country']}
    regions_current = {v['region'] for v in samples.values() if v['country']}
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
                    list(region_colors.columns)+[""]+list(countries_current), list(mat_reg[0])+[""]+geo_color])


    #write out colours
    data = data.transpose()
    color_info = pd.DataFrame(data)
    color_info = color_info.to_csv(sep='\t', index=False, header=False) #, header=['property', 'name', 'color'])
    with open(args.output_color, 'w') as outfile:
        outfile.write(color_info)

    #write out meta
    (pd.DataFrame.from_dict(data=samples, orient='index').to_csv(args.output_meta, sep='\t', header=True, index_label=False, index=False))

    # ...
    # update files: geo_info/geo_lat_long/properties_info/synonyms
    # ...
    #geo_regions
    (pd.DataFrame.from_dict(data=region_dict, orient='index').to_csv(args.geo_regions, sep='\t', header=False))

    #synonyms
    write_synonyms(synonyms, args.synonyms)

    #write out location dict
    write_out_loc_file = input("\nWould you like to save a file linking your locations and countries? [y] or n: ")
    if write_out_loc_file != 'n':
        geo_file_name = args.geo_regions.split("/")[-1]
        location_file = args.geo_regions.replace(geo_file_name, "location_country.tsv")
        (pd.DataFrame.from_dict(data=location_country, orient='index').to_csv(location_file, sep='\t', header=False))
        print("Location/country data written out to %s"%(location_file))

    #Print out lat-long file if user wants it - and if new coordinates have been added
    if len(coordinates) != orig_coordinate_length:
        print("\nNew location lat/long values were found.")
        write_out_coord = input("Would you like to write them out to a file? [y] or n: ")
        if write_out_coord != 'n':
            geo_file_name = args.geo_regions.split("/")[-1]
            location_file = args.geo_regions.replace(geo_file_name, "lat_longs.tsv")
            pd.DataFrame.from_dict(data=coordinates, orient='index').to_csv(location_file, sep='\t', header=False)
            print("Lat/long data written out to %s"%(location_file))


    #if include hosts/properties someday, do here

    #filtering and writing fasta would go here if included.
