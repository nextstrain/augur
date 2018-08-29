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

def id_standard_columns(orig_meta):
    #force user to keep required columns
    print('These are the available columns:')
    for i, j in enumerate(orig_meta.columns):
        print(i, j)
    print(i+1, 'none of them')
    take = []
    tc1 = orig_meta.columns[int(input('Which one of your columns is "strain" (matches sequence names)? '))]
    tc2 = orig_meta.columns[int(input('Which one of your columns is "accession"? '))]
    tc3 = orig_meta.columns[int(input('Which one of your columns is "date"? '))]
    tc4 = orig_meta.columns[int(input('Which one of your columns is "country"? '))]
    take.extend([tc1,tc2,tc3,tc4])

    return take

def id_keep_columns(orig_meta, columns_to_keep):
    print('\n\nRegion information will be generated within this script.')
    print('If you have region information in the data it will be in conflict with regions generated here.')
    print('One solution is to not take column "region".')
    print('These are the available columns not already taken:')
    for i, j in enumerate(orig_meta.columns):
        if j not in columns_to_keep:
            print('- ',j)
            
    for i, j in enumerate(orig_meta.columns):
        if j not in columns_to_keep:
            print("Take column \'{}\'?".format(j))
            c = input("[y] or n: ")
            
            if c.lower() != 'n':
                columns_to_keep.append(j)

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
def adjust_date(meta):

    drop_tsv = []
    print(meta.loc[10:35,'date'])
    print('Example date format: %m-%d-%Y')
    date_format = input("Enter date format: ")  
    minmax = [2000, 2100]  #set this to be something better 

    for s in meta.index:

        if len(meta.loc[s, 'date'])==4:
            meta.loc[s, 'date'] = meta.loc[s, 'date'] + '-XX-XX'
            if parse_date(meta.loc[s,'date'], '%Y-%m-%d', minmax)[0] == None:
                drop_tsv.append(s)

        elif len(meta.loc[s, 'date']) == 7: 
            if 'm' in date_format.split('Y')[0]:
                meta.loc[s, 'date'] = meta.loc[s, 'date'][3:] + '-' + meta.loc[s, 'date'][:2] + '-XX'                
            else:
                meta.loc[s, 'date'] = meta.loc[s, 'date'][:4] + '-' + meta.loc[s, 'date'][5:7] + '-XX'

            if parse_date(meta.loc[s,'date'], '%Y-%m-%d', minmax)[0] == None:
                drop_tsv.append(s)

        elif 'X' in meta.loc[s, 'date']:
            dec_date = parse_date(meta.loc[s, 'date'], date_format, minmax)
            if dec_date[0] == None:
                drop_tsv.append(s)
            else:            
                part1 = meta.loc[s,'date'].split(date_format.split('%')[1][1:2])[0]
                part2 = meta.loc[s,'date'].split(date_format.split('%')[1][1:2])[1]
                part3 = meta.loc[s,'date'].split(date_format.split('%')[1][1:2])[2]
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

                meta.loc[s, 'date'] = year + '-' + month + '-XX'

        else:    
            dec_date = parse_date(meta.loc[s,'date'], date_format, minmax)
            if dec_date == None:
                drop_tsv.append(s)
            else:
                meta.loc[s,'date'] = datetime.strptime(meta.loc[s,'date'], date_format).strftime("%Y-%m-%d")

    return(meta, drop_tsv)

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
def region_questionnaire(meta, region_dict, region_colors, a):
    
    if meta.loc[a,'country'] in region_dict:
        meta.loc[a, 'region'] = region_dict[meta.loc[a,'country']]                       
    else:
        print("Available regions are:")        
        for i, j in enumerate(region_colors.columns):
            print(i, j)

        print("Unknown location:", meta.loc[a,'country'])        
        reg_no = input("Type the region number to specify which region the country belongs to: ")    
        meta.loc[a, 'region'] = region_colors.columns[int(reg_no)]

        region_dict[meta.loc[a,'country']] = meta.loc[a, 'region']

    return(meta, region_dict)

def geo_color_generator(countries_current, region_counter, region_dict, region_colors):

    color = []

    for place in countries_current:
        html_col = []
        region = region_dict[place]
        area_diversity = 4 * region_counter[region]
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
        color.append(c_c)

    return(color)

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




def run(args):

    geolocator = Nominatim()

    #Read in meta-data provided by user
    orig_meta = pd.read_csv(args.input, sep="\t")

    #Get user to identify columns, select those to keep
    columns_to_keep = id_standard_columns(orig_meta)
    standard_columns = columns_to_keep[:] #copy to keep seperate

    #Pick other columns to keep
    columns_to_keep = id_keep_columns(orig_meta, columns_to_keep)

    ## properties section would go here - leaving out for now.
    ##### (this standardizes your host names, etc)
    #####

    non_standard_columns = [c for c in columns_to_keep if c not in standard_columns]

    new_meta = orig_meta.loc[:, columns_to_keep]
    new_meta['region'] = 0  


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
    new_meta, drop_tsv = adjust_date(new_meta)

    # ...            
    # check for duplicates
    # ...    
    # TODO: print some kind of message if removing some duplicates
    for i in range(len(new_meta.loc[:,'strain'])-1):
        if new_meta.loc[i, 'strain'] == new_meta.loc[i+1, 'strain']:
            drop_tsv.append(i+1)

    # ...
    # filter meta (drop duplicates and the ones with unparsable dates)
    # ...
    #TODO: what is 'unparsable dates'? Are those with no dates dropped?
    '''
    (note from sanda): this can have additional filtering criteria, just append to drop_tsv
    perhaps import and append indexes from dropped_strains.txt if you need them to be dropped as well
    '''   
    new_meta.drop(new_meta.index[drop_tsv], inplace=True)   



    ##Get region information - file assigns each country a region
    region_dict = {}
    with open(args.geo_regions, mode='r') as f:
        for line in f:
            vals = line.strip().split('\t')
            region_dict[vals[0]] = vals[1]
    
    region_colors = pd.read_csv(args.region_colors, sep='\t')
    region_counter = {}
    for r in region_colors.columns:
        region_counter[str(r)] = 0  

    ##Get country synonyms. We store in one format (more human friendly), and 'use' in another
    synonyms = get_synonyms(args.synonyms)

    # Get lat-longs - needs to be adapted to new format!!
    coor = pd.read_csv(args.geo_info, sep='\t') 
    #coor = read_lat_longs()  # should use this instead! But needs adapting

    #go through table...
    countries_current = [] #this keeps track for generating colours later...

    for a in new_meta.index:    
        wrong_name = ""
        #country = new_meta.loc[a, 'country']

        if type(new_meta.loc[a, 'country']) == str: 
           new_meta.loc[a, 'country'] = new_meta.loc[a, 'country'].lower().replace(' ','_')
        else:
            new_meta.loc[a, 'country'] = default_unknown

        if new_meta.loc[a, 'country'] in region_dict:
            new_meta.loc[a, 'region'] = region_dict[new_meta.loc[a, 'country']]   

        elif new_meta.loc[a, 'country'] in indicators_for_unknown:        
            new_meta.loc[a, 'country'] = default_unknown
            new_meta.loc[a,'region'] = default_unknown
                 
        elif new_meta.loc[a, 'country'] in synonyms:
            new_meta.loc[a, 'country'] = synonyms[new_meta.loc[a, 'country']]
            new_meta.loc[a, 'region'] = region_dict[new_meta.loc[a, 'country']]   

        else:
            wrong_name = new_meta.loc[a, 'country']
            location = geolocator.geocode((wrong_name).replace('_',' '), language='en')
            if str(location) == 'None':
                print("The location as formatted couldn't be found in GeoPy")
                answer = 'n'
            else:
                new_loc = (location.address).lower().replace(' ','_') 
                print("the unknown location is: ", wrong_name)
                print("GeoPy suggests using:", new_loc)
                answer = input("Is this the right country name? y or n: ")                
            
            if answer.lower() == 'y':
                new_meta.loc[a, 'country'] = new_loc
                if new_loc != list(coor.location[:]):
                    lon = location.longitude
                    lat = location.latitude
                    coor.loc[len(coor)] = [new_loc, 'XX', lat, lon]
                
                new_meta, region_dict = region_questionnaire(new_meta, region_dict, region_colors,  a)

            else:
                # Here we let the user try to type something else to find in synonyms or geoPy
                redo = True
                while redo == True:
                    print(wrong_name)
                    new_loc = (input("Type the correct country name or 'n' if unknown: ")).lower().replace(' ','_')
                    if new_loc in ['n', 'na']:
                        indicators_for_unknown.append(wrong_name)
                        new_meta.loc[a, 'country'] = default_unknown
                        new_meta.loc[a, 'region'] = default_unknown

                    else:
                        new_loc = new_loc.lower().replace(' ','_')
                        if new_loc in synonyms.keys():
                            new_loc = synonyms[new_loc]

                        new_meta.loc[a, 'country'] = new_loc

                        if new_loc in region_dict: #what was typed, or synonym of what was typed, is in geo
                            new_meta.loc[a, 'region'] = region_dict[new_loc]  
                            redo = False
                        else:
                            if new_loc != 'n' and new_loc in list(coor.location[:]):
                                new_meta, region_dict = region_questionnaire(new_meta, region_dict, region_colors,  a)
                                redo = False
                        
                            elif new_loc != 'n' and new_loc not in list(coor.location[:]):                  
                                location = geolocator.geocode((new_loc).replace('_',' '), language='en')
                                print('would you like to: \n1. Use the coordinates for this location?', location.address)
                                print('2. Enter coordinates manually in the file after this run')
                                print('3. Type the country again')
                                loc_ok = input('1, 2, or 3: ')
                                if loc_ok == '1':
                                    lon = location.longitude
                                    lat = location.latitude                    
                                    coor.loc[len(coor)] = [new_loc, 'XX', lat, lon]
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
                                        new_meta.loc[a, 'country'] = (location.address).lower().replace(' ','_') 
                                        new_loc = new_meta.loc[a,'country']
                                    new_meta, region_dict = region_questionnaire(new_meta, region_dict, region_colors, a)


            print(wrong_name)
            s_add = input("Was this country misspelled (should it be added to synonyms)? y or n: ")
            if s_add.lower() == 'y':
                synonyms[wrong_name] = new_loc  

        if new_meta.loc[a,'country'] != default_unknown:        
            if new_meta.loc[a,'country'] not in countries_current:
                countries_current.append(new_meta.loc[a,'country'])
                region_counter[new_meta.loc[a, 'region']]+=1

    ## Properties things will go here eventually (host, etc)
    ####
    ####

    # ...
    # generate region/country colors
    # ...
    geo_color = geo_color_generator(countries_current, region_counter, region_dict, region_colors)


    
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
