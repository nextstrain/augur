"""
all of this information can one day make it to the DB I would think
It lives in a seperate file simply to make flu.prepare.py less cluttered
"""

segments = ["pb2", "pb1", "pa", "ha", "np", "na", "mp", "ns"]

# regions is list of tuples (region, acronym)
# acronym = "" means ignore for frequency calcs
regions = [
    ('africa',          ""),
    ('south_asia',      "AS"),
    ('europe',          "EU"),
    ('china',           "AS"),
    ('north_america',   "NA"),
    ('south_america',   ""),
    ('japan_korea',     "AS"),
    ('oceania',         "OC"),
    ('southeast_asia',  "AS"),
]

outliers = {
    'h3n2':["A/Sari/388/2006", "A/SaoPaulo/36178/2015", "A/Pennsylvania/40/2010", "A/Pennsylvania/14/2010",
            "A/Pennsylvania/09/2011", "A/OSAKA/31/2005", "A/Ohio/34/2012", "A/Kenya/170/2011", "A/Kenya/168/2011",
            "A/Indiana/21/2013", "A/Indiana/13/2012", "A/Indiana/6/2013", "A/Indiana/17/2013", "A/Indiana/11/2013", "A/Indiana/08/2012", "A/Indiana/06/2013",
            "A/India/6352/2012", "A/HuNan/01/2014", "A/Helsinki/942/2013", "A/Guam/AF2771/2011", "A/Chile/8266/2003",
            "A/Busan/15453/2009", "A/Nepal/142/2011", "A/Kenya/155/2011", "A/Guam/AF2771/2011", "A/Michigan/82/2016",
            "A/Ohio/27/2016", "A/Ohio/28/2016", "A/Michigan/83/2016", "A/Michigan/84/2016", "A/Jiangsu-Tianning/1707/2013",
            "A/HuNan/1/2014", "A/Iran/227/2014", "A/Iran/234/2014", "A/Iran/140/2014", "A/Jiangsu-Chongchuan/1830/2014",
            "A/Chile/8266/2003", "A/Louisiana/4/2003", "A/Lousiana/4/2003", "A/OSAKA/31/2005",
            "A/Sari/388/2006", "A/HongKong/HK1/2008", "A/HongKong/HK1MA21-1/2008","A/HongKong/HK1MA21-4/2008", "A/HongKong/HK1MA21-2/2008",
            "A/HongKong/HK1MA21-3/2008", "A/HongKong/HK2/2008", "A/HongKong/HK2MA21-1/2008",
            "A/HongKong/HK2MA21-2/2008", "A/HongKong/HK2MA21-3/2008", "A/HongKong/HK4/2008",
            "A/HongKong/HK5/2008", "A/HongKong/HK5MA21-1/2008", "A/HongKong/HK5MA21-3/2008",
            "A/HongKong/HK6/2008", "A/HongKong/HK6MA21-2/2008", "A/HongKong/HK6MA21-3/2008",
            "A/HongKong/HK7/2008", "A/HongKong/HK8/2008", "A/HongKong/HK8MA21-1/2008",
            "A/HongKong/HK8MA21-2/2008", "A/HongKong/HK8MA21-3/2008", "A/HongKong/HK8MA21-4/2008",
            "A/HongKong/HK9/2008", "A/HongKong/HK9MA21-1/2008", "A/HongKong/HK9MA21-2/2008",
            "A/HongKong/HK9MA21-3/2008", "A/HongKong/HK10/2008", "A/HongKong/HK10MA21-1/2008",
            "A/HongKong/HK10MA21-2/2008", "A/HongKong/HK10MA21-3/2008", "A/HongKong/HK10MA21-4/2008",
            "A/HongKong/HK11MA21-1/2008", "A/HongKong/HK11MA21-3/2008", "A/HongKong/HK11MA21-4/2008",
            "A/HongKong/HK12/2008", "A/HongKong/HK12MA21-2/2008", "A/HongKong/HKMA12/2008",
            "A/HongKong/HKMA12A/2008", "A/HongKong/HKMA12B/2008", "A/HongKong/HKMA12D/2008",
            "A/HongKong/HKMA12E/2008", "A/HongKong/HKMA20B/2008", "A/HongKong/HKMA20E/2008", "A/Kansas/13/2009",
            "A/Busan/15453/2009", "A/Pennsylvania/14/2010", "A/Pennsylvania/40/2010", "A/Guam/AF2771/2011",
            "A/Indiana/8/2011", "A/Kenya/155/2011", "A/Kenya/168/2011", "A/Kenya/170/2011", "A/Nepal/142/2011",
            "A/Pennsylvania/09/2011", "A/Pennsylvania/9/2011", "A/Quebec/167936/2011", "A/Quebec/170658/2011",
            "A/India/6352/2012", "A/Indiana/08/2012", "A/Indiana/13/2012", "A/Ohio/34/2012",
            "A/Helsinki/942/2013", "A/Indiana/06/2013", "A/Indiana/11/2013", "A/Indiana/21/2013",
            "A/Jiangsu-Tianning/1707/2013", "A/HuNan/01/2014", "A/Jiangsu-Chongchuan/1830/2014",
            "A/Jiangsu-Chongchuan/12179/2014", "A/Ohio/2/2014", "A/Ohio/4319/2014", "A/SaoPaulo/3-34891/2014",
            "A/Wisconsin/24/2014", "A/NewJersey/53/2015", "A/SaoPaulo/36178/2015", "A/SaoPaulo/61282/2015",
            "A/SaoPaulo/65194/2015", "A/Michigan/39/2015", "A/Sydney/53/2015", "A/Michigan/82/2016",
            "A/Michigan/83/2016", "A/Michigan/84/2016", "A/Michigan/87/2016", "A/Michigan/89/2016",
            "A/Michigan/90/2016", "A/Michigan/91/2016", "A/Michigan/93/2016", "A/Michigan/94/2016",
            "A/Michigan/95/2016", "A/Michigan/96/2016", "A/Ohio/27/2016", "A/Ohio/28/2016", "A/Ohio/32/2016",
            "A/Ohio/33/2016", "A/Ohio/35/2016", "A/Zhejiang-Wuxin/1300/2016", "A/Nanjing/1/2010", "A/StPetersburg/5/2009",
            "A/Cambodia/NHRCC00001/2009","A/Cambodia/NHRCC00002/2009", "A/Cambodia/NHRCC00003/2009",
            "A/Cambodia/NHRCC00006/2009", "A/Iran/407/2014"],
    'h1n1pdm': ["A/Kenya/264/2012", "A/Iowa/39/2015", "A/Asturias/RR6898/2010", "A/Wisconsin/28/2011",
                "A/Brest/1161/2014", "A/Tomsk/273-MA1/2010", "A/Minnesota/46/2015", "A/Poland/16/2013",
                "A/Hungary/02/2013", "A/Hungary/16/2013", "A/California/07/2009NYMC-X18113/198",
                "A/Christchurch/16/2010NIB-74xp13/202", "A/Bari/166/2016", "A/Bari/167/2016", "A/Dakar/3/2014",
                "A/Arkansas/15/2013", "A/Wisconsin/87/2005", 'A/Norway/1929/2014', 'A/Ohio/9/2015'],
    'vic':["B/Bangkok/SI17/2012", "B/Bangkok/SI58/2012", "B/Kol/2024/2008", "B/Kolkata/2024/2008"],
    "yam":[]
}

reference_maps = {
    "h3n2": {
        "ha": {
            "path": "metadata/h3n2_ha_outgroup.gb",
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
            "use": True,
            "genes": ["HA1", "HA2"]
        }
    },
    "yam": {
        "ha": {
            "path": "metadata/yam_ha_outgroup.gb",
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
            "use": True,
            "genes": ["HA"]
        },
        "na": {
            "path": "metadata/yam_na_outgroup.gb",
            "use": True,
            "genes": ["NA", "NB"]
        }
    },
    "vic": {
        "ha": {
            "path": "metadata/vic_ha_outgroup.gb",
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
            "use": True,
            "genes": ["HA"]
        }
    },
    "h1n1pdm": {
        "ha": {
            "path": "metadata/h1n1pdm_ha_outgroup.gb",
            "metadata": {
                'strain': "A/Swine/Indiana/P12439/00",
                'isolate_id': "AF455680",
                'date': "unknown",
                'region': "north america",
                'country': "USA",
                "city": "unknown",
                "passage": "unknown",
                'lab': "unknown",
                'age': "unknown",
                'gender': "unknown"
            },
            "use": False,
            "genes": ["HA"]
        }
    }
}

## lots of the references share data
reference_maps["yam"]["na"]["metadata"] = reference_maps["yam"]["ha"]["metadata"]

reference_viruses = {
    'h3n2': ['A/Wisconsin/67/2005', 'A/Brisbane/10/2007',  'A/Perth/16/2009', 'A/Victoria/361/2011','A/Texas/50/2012', 'A/Switzerland/9715293/2013', 'A/HongKong/4801/2014', 'A/Alaska/232/2015'],
    'h1n1pdm':[],
    'vic':[],
    'yam':[]
}
