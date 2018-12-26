"""
all of this information can one day make it to the DB I would think
It lives in a seperate file simply to make flu.prepare.py less cluttered
"""

segments = ["pb2", "pb1", "pa", "ha", "np", "na", "mp", "ns"]

# regions is list of tuples (region, acronym, popsize in billions)
# acronym = "" means ignore for frequency calcs
regions = [
    ('africa',            "",   1.02),
    ('europe',            "EU", 0.74),
    ('north_america',     "NA", 0.54),
    ('china',             "AS", 1.36),
    ('south_asia',        "AS", 1.45),
    ('japan_korea',       "AS", 0.20),
    ('oceania',           "OC", 0.04),
    ('south_america',     "SA", 0.41),
    ('southeast_asia',    "AS", 0.62),
    ('west_asia',         "AS", 0.75)
]

reference_maps = {
    "h3n2": {
        "ha": {
            "path": "metadata/h3n2_outgroup.gb",
            "metadata": {
                'strain': "A/Beijing/32/1992",
                'isolate_id': "CY113677",
                'date': "1992-XX-XX",
                'region': "china",
                'country': "China",
                "city": "Beijing",
                "passage": "unknown",
                'lab': "unknown",
                'age': "unknown",
                'gender': "unknown"
            },
            "include": 0,
            "genes": ["SigPep", "HA1", "HA2"]
        },
        "na": {
            "path": "metadata/h3n2_na_outgroup.gb",
            "metadata": {
                'strain': "A/Beijing/32/1992",
                'isolate_id': "CY113677",
                'date': "1992-XX-XX",
                'region': "china",
                'country': "China",
                "city": "Beijing",
                "passage": "unknown",
                'lab': "unknown",
                'age': "unknown",
                'gender': "unknown"
            },
            "include": 0,
            "genes": ["NA"]
        }
    },
    "yam": {
        "ha": {
            "path": "metadata/yam_outgroup.gb",
            "metadata": {
                'strain': "B/Singapore/11/1994",
                'isolate_id': "CY019707",
                'date': "1994-05-10",
                'region': "southeast_asia",
                'country': "Singapore",
                "city": "Singapore",
                "passage": "unknown",
                'lab': "unknown",
                'age': "unknown",
                'gender': "M"
            },
            "include": 0,
            "genes": ["SigPep", "HA1", "HA2"]
        },
        "na": {
            "path": "metadata/yam_na_outgroup.gb",
            "metadata": {
                'strain': "B/Singapore/11/1994",
                'isolate_id': "CY019707",
                'date': "1994-05-10",
                'region': "southeast_asia",
                'country': "Singapore",
                "city": "Singapore",
                "passage": "unknown",
                'lab': "unknown",
                'age': "unknown",
                'gender': "M"
            },
            "include": 0,
            "genes": ["NA", "NB"]
        }
    },
    "vic": {
        "ha": {
            "path": "metadata/vic_outgroup.gb",
            "metadata": {
                'strain': "B/Hong Kong/02/1993",
                'isolate_id': "CY018813",
                'date': "1993-02-15",
                'region': "china",
                'country': "Hong Kong",
                "city": "Hong Kong",
                "passage": "4",
                'lab': "unknown",
                'age': "unknown",
                'gender': "unknown"
            },
            "include": 0,
            "genes": ["SigPep", "HA1", "HA2"]
        },
        "na": {
            "path": "metadata/vic_na_outgroup.gb",
            "metadata": {
                'strain': "B/Hong Kong/02/1993",
                'isolate_id': "CY018813",
                'date': "1993-02-15",
                'region': "china",
                'country': "Hong Kong",
                "city": "Hong Kong",
                "passage": "4",
                'lab': "unknown",
                'age': "unknown",
                'gender': "unknown"
            },
            "include": 0,
            "genes": ["NB", "NA"]
        }
    },
    "h1n1pdm": {
        "ha": {
            "path": "metadata/h1n1pdm_ha_outgroup.gb",
            "metadata": {
                'strain': "A/California/07/2009",
                'isolate_id': "CY121680",
                'date': "2009-04-09",
                'region': "north america",
                'country': "USA",
                "city": "unknown",
                "passage": "egg",
                'lab': "unknown",
                'age': "54",
                'gender': "male"
            },
            "include": 0,
            "genes": ["SigPep", "HA1", "HA2"]
        },
        "na": {
            "path": "metadata/h1n1pdm_na_outgroup.gb",
            "metadata": {
                'strain': "A/California/07/2009",
                'isolate_id': "CY121682",
                'date': "2009-04-09",
                'region': "north america",
                'country': "USA",
                "city": "unknown",
                "passage": "egg",
                'lab': "unknown",
                'age': "54",
                'gender': "male"
            },
            "include": 0,
            "genes": ["NA"]
        }
    }
}

## lots of the references share data
reference_maps["yam"]["na"]["metadata"] = reference_maps["yam"]["ha"]["metadata"]

# these will be included in narrower builds, for example 'A/Michigan/15/2014' will appear in a 2015-2017
# place vaccine strains here - this ensures they'll be sampled
reference_viruses = {
    'h3n2':[
        'A/Wisconsin/67/2005', 'A/Brisbane/10/2007', 'A/Perth/16/2009', 'A/Victoria/361/2011',
        'A/Texas/50/2012', 'A/Southauckland/7/2012' 'A/Singapore/22/2012', 'A/Christchurch/516/2012',
        'A/Hawaii/22/2012', 'A/Switzerland/9715293/2013', 'A/Victoria/505/2013', 'A/Victoria/506/2013',
        'A/Almaty/2958/2013', 'A/HongKong/4801/2014', 'A/Hawaii/47/2014', 'A/Tasmania/11/2014',
        'A/Newcastle/22/2014', 'A/Sidney/7/2014', 'A/Netherlands/524/2014', 'A/Norway/466/2014',
        'A/Michigan/15/2014', 'A/NewJersey/26/2014', 'A/HongKong/7127/2014', 'A/SouthAustralia/55/2014',
        'A/NewCaledonia/71/2014', 'A/Stockholm/6/2014', 'A/Fiji/2/2015', 'A/Wisconsin/20/2015',
        'A/Wisconsin/66/2015', 'A/Brisbane/45/2015', 'A/Ontario/RV2414/2015', 'A/Nebraska/19/2015',
        'A/Montana/28/2015', 'A/Alaska/232/2015', 'A/Alaska/240/2015', 'A/Singapore/Infimh-16-0019/2016',
        'A/Texas/88/2016', 'A/Idaho/33/2016', 'A/NorthCarolina/4/2016', 'A/Delaware/32/2016',
        'A/Nevada/22/2016', 'A/Brisbane/321/2016', 'A/Newcastle/30/2016', 'A/Norway/3806/2016',
        'A/Canberra/7/2016', 'A/Washington/106/2016', 'A/Brisbane/318/2016', 'A/Oregon/19/2016',
        'A/Antananarivo/1067/2016', 'A/Florida/23/2017', 'A/Washington/16/2017', 'A/NorthCarolina/4/2017',
        'A/SouthCarolina/4/2017', 'A/Texas/71/2017', 'A/Wisconsin/327/2017', 'A/Texas/68/2017',
        'A/Brisbane/29/2017', 'A/Nebraska/2/2017',
        'A/Texas/50/2012-egg', 'A/Switzerland/9715293/2013-egg', 'A/HongKong/7127/2014-egg',
        'A/NewCaledonia/71/2014-egg', 'A/HongKong/4801/2014-egg', 'A/Saitama/103/2014-egg',
        'A/Alaska/232/2015-egg', 'A/Singapore/Infimh-16-0019/2016-egg', 'A/Norway/3806/2016-egg',
        'A/Idaho/33/2016-egg', 'A/Antananarivo/1067/2016-egg', 'A/HongKong/50/2016-egg',
        'A/Norway/4465/2016-egg', 'A/Wisconsin/19/2017-egg', 'A/Greece/4/2017-egg',
        'A/HongKong/2286/2017-egg', 'A/Victoria/653/2017-egg', 'A/Norway/4465/2017', 'A/Norway/4465/2017-egg',
        'A/Switzerland/8060/2017-egg', 'A/Bretagne/1565/2017-egg', 'A/Norway/4436/2016', 'A/Switzerland/8060/2017'
    ],
    'h1n1pdm':[
        'A/California/7/2009', 'A/Bangladesh/2021/2012', 'A/Victoria/367/2012', 'A/Brisbane/96/2012',
        'A/SouthAfrica/3626/2013', 'A/Brisbane/28/2013', 'A/Bolivia/559/2013', 'A/Florida/62/2014',
        'A/Newcaledonia/58/2014', 'A/Tasmania/24/2014', 'A/Florida/62/2014', 'A/Michigan/45/2015',
        'A/Perth/103/2015', 'A/Brisbane/181/2015', 'A/Singapore/GP1908/2015', 'A/Iowa/53/2015',
        'A/StPetersburg/61/2015', 'A/Minnesota/32/2015', 'A/Bangladesh/3002/2015', 'A/Israel/Q-504/2015',
        'A/Indiana/21/2016', 'A/Victoria/503/2016', 'A/Fiji/3/2016', 'A/Panama/318595/2016',
        'A/Montana/50/2016', 'A/Texas/157/2016', 'A/Michigan/272/2017', 'A/Perth/10/2017',
        'A/Arizona/33/2017', 'A/Brisbane/37/2017', 'A/Oman/3192/2017', 'A/AbuDhabi/18/2017',
        'A/India/3405/2017', 'A/India/322/2017', 'A/Navarra/1985/2017', 'A/Kazakhstan/581/2018'
    ],
    'vic':[
        'B/Shangdong/7/1997', 'B/HongKong/330/2001', 'B/Malaysia/2506/2004', 'B/Ohio/1/2005',
        'B/Brisbane/60/2008', 'B/Utah/8/2012', 'B/Montana/5/2012', 'B/Texas/2/2013', 'B/Florida/33/2014',
        'B/Indiana/25/2015', 'B/Florida/78/2015', 'B/Iceland/56/2015', 'B/Missouri/9/2016',
        'B/Ireland/3154/2016', 'B/StPetersburg/293/2016', 'B/Maryland/15/2016', 'B/Colorado/6/2017',
        'B/Kazakhstan/128/2017', 'B/HongKong/286/2017', 'B/Bangladesh/3004/2017',
        'B/Vologda/CRIE-323/2016', 'B/NewJersey/4/2017', 'B/Indiana/25/2015', 'B/Singapore/TT1243/2015', 'B/Brisbane/312/2015'
    ],
    'yam':[
        'B/Beijing/184/1993', 'B/Sichuan/379/1999', 'B/Shanghai/361/2002', 'B/Florida/4/2006',
        'B/Wisconsin/1/2010', 'B/Massachusetts/2/2012', 'B/Phuket/3073/2013', 'B/Utah/9/2014',
        'B/Brisbane/9/2014', 'B/Guangdong-Liwan/1133/2014', 'B/California/12/2015', 'B/Arizona/10/2015',
        'B/Sapporo/2/2015', 'B/Hyogo/3210/2015', 'B/NewHampshire/1/2016', 'B/Texas/81/2016'
    ]
}

vaccine_choices = {
    "h3n2": {
        'A/Sydney/5/1997': "1997-09-25",
        'A/Moscow/10/1999': "1999-09-25",
        'A/Fujian/411/2002': "2003-09-25",
        'A/California/7/2004': "2005-02-21",
        'A/Wisconsin/67/2005': "2006-02-21",
        'A/Brisbane/10/2007': "2007-09-25",
        'A/Perth/16/2009': "2009-09-25",
        'A/Victoria/361/2011': "2012-02-21",
        'A/Texas/50/2012': "2013-09-25",
        'A/Switzerland/9715293/2013': "2014-09-25",
        'A/HongKong/4801/2014': "2015-09-24",
        'A/Singapore/Infimh-16-0019/2016': "2017-09-28",
        'A/Switzerland/8060/2017': "2018-09-27"
    },
    "vic": {
        'B/Shangdong/7/1997': "1999-09-25",
        'B/HongKong/330/2001': "2002-09-25",
        'B/Malaysia/2506/2004': "2006-09-25",
        'B/Brisbane/60/2008': "2009-09-25",
        'B/Colorado/6/2017': "2018-02-22"
    },
    "yam": {
        'B/Beijing/184/1993': "1998-11-01",
        'B/Sichuan/379/1999': "2001-09-25",
        'B/Shanghai/361/2002': "2004-09-25",
        'B/Florida/4/2006': "2008-09-25",
        'B/Wisconsin/1/2010': "2012-02-25",
        'B/Massachusetts/2/2012': "2013-02-25",
        'B/Phuket/3073/2013': "2014-09-25"
    },
    "h1n1pdm": {
        'A/California/7/2009': "2009-09-25",
        'A/Michigan/45/2015': "2016-09-29"
    }
}

# Local Branching Index (LBI) params
LBI_params = {
    '2y': {"tau": 0.3, "time_window": 0.5},
    '3y': {"tau": 0.4, "time_window": 0.6},
    '6y': {"tau": 0.25, "time_window": 0.75},
    '12y': {"tau": 0.25, "time_window": 0.75}
}

# Frequency Params
frequency_params = {
    '2y': {"dfreq_dn": 6},
    '3y': {"dfreq_dn": 6},
    '6y': {"dfreq_dn": 6},
    '12y': {"dfreq_dn": 6}
}

# Map resolution to pivot spacing
resolution_to_pivot_spacing = {
    "2y": 1. / 12.,
    "3y": 1. / 12.,
    "6y": 2. / 12.,
    "12y": 3. / 12.
}

# Map lineages to specific HA masks
lineage_to_epitope_mask = {
    "h3n2": "wolf",
    "h1n1pdm": "canton"
}

# Map lineages to specific NA masks
lineage_to_na_epitope_mask = {
    "h3n2": "bhatt"
}

lineage_to_glyc_mask = {
    "h3n2": "ha1_h3n2",
    "h1n1pdm": "ha1_globular_head_h1n1pdm"
}

clade_designations = {
    "h3n2": {
        '3b':   [('HA1',48,'T'), ('HA2',158,'N'), ('HA1',198,'S'), ('HA1',312,'S'), ('HA1',223,'I'), ('HA1',145,'S')],
        '3c':   [('HA1',45,'N'), ('HA1',48,'I'), ('nuc',456,'T'), ('HA1',198,'S'), ('HA1',312,'S'),
                 ('HA1',223,'I')],
        '3c3':  [('nuc',285,'T'), ('nuc',430,'G'), ('nuc',472,'G'), ('nuc',1296,'A')],
        '3c3.A':[('HA1',128,'A'), ('HA1',142,'G'), ('HA1',159,'S')],
        '3c3.B':[('HA1',83,'R'), ('HA1',261,'Q'), ('HA1',62,'K'), ('HA1',122,'D')],
        '3c2':  [('HA2',160,'N'), ('HA1',145,'S'), ('nuc',693,'A'), ('nuc',1518,'G')],
        '3c2.A':[('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'), ('HA1',311,'H'), ('nuc',1491,'A'),
                 ('nuc',234,'A')],
        'A1':   [('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'), ('HA1',311,'H'), ('nuc',1491,'A'),
                 ('HA1',171,'K'), ('HA2',77,'V'), ('HA2',155,'E'), ('HA2',160,'N')],
        'A2':   [('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'), ('HA1',311,'H'), ('nuc',1491,'A'),
                 ('nuc',234,'A'), ('HA1',131,'K'), ('HA1',142,'K'), ('HA1',261,'Q')],
        'A2/re':[('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'), ('HA1',311,'H'), ('nuc',1491,'A'),
                 ('nuc',234,'A'), ('HA1',131,'K'), ('HA1',142,'K'), ('HA1',261,'Q'), ('nuc',1689,'T')],
        'A3':   [('HA1',159,'Y'), ('HA1',225,'D'), ('HA1',311,'H'), ('nuc',1491,'A'), ('nuc',234,'A'),
                 ('HA1',121,'K'), ('HA1',144,'K')],
        'A4':   [('HA1',159,'Y'), ('HA1',225,'D'), ('HA1',311,'H'), ('nuc',1491,'A'), ('nuc',234,'A'),
                 ('HA1',53,'N'), ('HA1',144,'R'), ('HA1',171,'K'), ('HA1',192,'T'), ('HA1',197,'H')],
        'A1a':  [('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'), ('HA1',311,'H'), ('nuc',1491,'A'),
                 ('nuc',234,'G'), ('HA2',150,'E'), ('nuc',114,'T')],
        'A1b':  [('HA1',159,'Y'), ('HA1',225,'D'), ('nuc',1491,'A'), ('nuc',234,'G'), ('HA1',92,'R'),
                 ('HA1',311,'Q'), ('nuc',538,'C')],
        'A1b/135N': [('HA1',159,'Y'), ('HA1',225,'D'), ('nuc',1491,'A'), ('nuc', 234,'G'), ('HA1',92,'R'),
                     ('HA1',311,'Q'), ('nuc',538,'C'), ('nuc',81,'G'), ('nuc',453,'T')],
        'A1b/135K': [('HA1',159,'Y'), ('HA1',225,'D'), ('nuc',1491,'A'), ('nuc', 234, 'G'), ('HA1',92,'R'),
                     ('HA1',311,'Q'), ('nuc',538,'C'), ('nuc',233,'G'), ('nuc',472,'G'), ('nuc',452,'A')]
    },
    "h1n1pdm": {
        '1':        [('HA1',125,'N'), ('HA1',134,'A'), ('HA1',183,'S'), ('HA1',31,'N'), ('HA1',216,'I')],
        '2':        [('HA1',125,'N'), ('HA1',256,'A'), ('HA1',134,'A'), ('HA1',183,'S'), ('HA1',31,'D')],
        '3':        [('HA1',134,'T'), ('HA1',183,'P')],
        '4':        [('HA1',125,'D'), ('HA1',134,'A'), ('HA1',183,'S')],
        '5':        [('HA1',87,'N'), ('HA1',205,'K'), ('HA1',216,'V'), ('HA1',149,'L')],
        '6':        [('HA1',185,'T'), ('HA1',97,'N'), ('HA1',197,'A')],
        '6c':       [('HA1',234,'I'), ('HA1',256,'A'), ('HA1',97,'N'), ('HA1',197,'A'), ('HA1',283,'E')],
        '6b':       [('HA1',163,'Q'), ('HA1',256,'T'), ('HA1',197,'A'), ('HA1',283,'E')],
        '7':        [('HA1',143,'G'), ('HA1',97,'D'), ('HA1',197,'T')],
        '8':        [('HA1',186,'T'), ('HA1',272,'A')],
        '6b1':      [('HA1',163,'Q'), ('HA1',256,'T'), ('HA1',197,'A'), ('HA1',283,'E'), ('SigPep',13, 'T'),
                     ('HA1',84,'N'), ('HA1',162,'N')],
        '6b2':      [('HA1',163,'Q'), ('HA1',256,'T'), ('HA1',197,'A'), ('HA1',283,'E'), ('HA2',164,'G'),
                     ('HA1',152,'T'), ('HA2',174,'E')],
        '6b1.A':    [('HA1',84,'N'), ('HA1',162,'N'), ('HA1',74,'R'), ('HA1',164,'T'), ('HA1',295,'V')],
        '6b1.A1':   [('HA1',84,'N'), ('HA1',162,'N'), ('HA1',74,'R'), ('HA1',164,'T'), ('HA1',295,'V'),
                     ('nuc',364,'C'), ('nuc',990,'G'), ('HA1',183,'P')],
        '6b1.A2':   [('HA1',84,'N'), ('HA1',162,'N'), ('HA1',74,'R'), ('HA1',164,'T'), ('HA1',295,'V'),
                     ('nuc',1344,'C'), ('nuc',1467,'C'), ('HA1',183,'P')],
        '6b1.A3':   [('HA1',84,'N'), ('HA1',162,'N'), ('HA1',74,'R'), ('HA1',164,'T'), ('HA1',295,'V'),
                     ('nuc',351,'A'), ('nuc',492,'T'), ('nuc',1206,'A'), ('HA1',120,'A')],
        '6b1.A4':   [('HA1',84,'N'), ('HA1',162,'N'), ('HA1',74,'R'), ('HA1',164,'T'), ('HA1',295,'V'),
                     ('nuc',1344,'T'), ('nuc',990,'A'), ('HA1',183,'P')],
        '6b1.A5':   [('HA1',84,'N'), ('HA1',162,'N'), ('HA1',74,'R'), ('HA1',164,'T'), ('HA1',295,'V'),
                     ('HA1',161,'I'), ('HA2',77,'M')],
    },
    "vic": {
        '1A':       [('nuc',206,'G'), ('nuc',644,'C'), ('nuc',1340,'T'), ('nuc',1821,'T'), ('HA1',165,'K'),
                     ('HA1',172,'P')],
        '1B':       [('nuc',1034,'G'), ('nuc',1172,'G'), ('HA1',165,'K'), ('HA1',172,'P')],
        'V1A':      [('HA1',75,'K'), ('HA1',58,'L'), ('HA1',165,'K'), ('HA1',129,'D'), ('HA1',117,'V')],
        'V1A.1':    [('HA1',75,'K'), ('HA1',58,'L'), ('HA1',165,'K'), ('HA1',117,'V'), ('HA1',180,'V'),
                     ('HA2',152,'K'), ('nuc',594,'G')],
        'V1A/165N': [('HA1',75,'K'), ('HA1',58,'L'), ('HA1',129,'D'), ('HA1',117,'V'), ('HA1',221,'I'),
                     ('HA1',165,'N')]
    },
    "yam": {
        '2':    [('HA1',48,'K'), ('HA1',108,'A'), ('HA1',150,'S')],
        '3':    [('HA1',48,'R'), ('HA1',108,'P'), ('HA1',150,'I')],
        '3a':   [('HA1',37,'A'), ('HA1',298,'E'), ('HA1',48,'R'), ('HA1',105,'P'), ('HA1',150,'I')],
        '172Q': [('HA1',48,'R'), ('HA1',108,'P'), ('HA1',150,'I'), ('HA1',116,'K'), ('nuc',848,'G'),
                 ('HA1',172,'Q'), ('HA1',298,'E'), ('HA1',312,'K')]
    }
}
