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
    ('south_america',     "",   0.41),
    ('southeast_asia',    "AS", 0.62),
    ('west_asia',         "AS", 0.75)
]

outliers = {
    'h3n2':[
        "A/Chile/8266/2003", "A/Louisiana/4/2003", "A/Lousiana/4/2003", "A/India/D0512577/2005",
        "A/Kunming/1-Va10/2005", "A/Ontario/RV1273/2005", "A/OSAKA/31/2005", "A/Sari/388/2006",
        "A/Santiago/6764/2006", "A/Ontario/1252/2007", "A/HongKong/HK1/2008", "A/HongKong/HK1MA21-1/2008",
        "A/HongKong/HK1MA21-2/2008", "A/HongKong/HK1MA21-3/2008", "A/HongKong/HK1MA21-4/2008",
        "A/HongKong/HK2/2008", "A/HongKong/HK2MA21-1/2008", "A/HongKong/HK2MA21-2/2008",
        "A/HongKong/HK2MA21-3/2008", "A/HongKong/HK4/2008", "A/HongKong/HK4MA21-3/2008",
        "A/HongKong/HK5/2008", "A/HongKong/HK5MA21-1/2008", "A/HongKong/HK5MA21-3/2008",
        "A/HongKong/HK6/2008", "A/HongKong/HK6MA21-1/2008", "A/HongKong/HK6MA21-2/2008",
        "A/HongKong/HK6MA21-3/2008", "A/HongKong/HK7/2008", "A/HongKong/HK8/2008",
        "A/HongKong/HK8MA21-1/2008", "A/HongKong/HK8MA21-2/2008", "A/HongKong/HK8MA21-3/2008",
        "A/HongKong/HK8MA21-4/2008", "A/HongKong/HK9/2008", "A/HongKong/HK9MA21-1/2008",
        "A/HongKong/HK9MA21-2/2008", "A/HongKong/HK9MA21-3/2008", "A/HongKong/HK10/2008",
        "A/HongKong/HK10MA21-1/2008", "A/HongKong/HK10MA21-2/2008", "A/HongKong/HK10MA21-3/2008",
        "A/HongKong/HK10MA21-4/2008", "A/HongKong/HK11/2008", "A/HongKong/HK11MA21-1/2008",
        "A/HongKong/HK11MA21-2/2008", "A/HongKong/HK11MA21-3/2008", "A/HongKong/HK11MA21-4/2008",
        "A/HongKong/HK12/2008", "A/HongKong/HK12MA21-2/2008", "A/HongKong/HKMA12/2008",
        "A/HongKong/HKMA12A/2008", "A/HongKong/HKMA12B/2008", "A/HongKong/HKMA12C/2008",
        "A/HongKong/HKMA12D/2008", "A/HongKong/HKMA12E/2008", "A/HongKong/HKMA20A/2008",
        "A/HongKong/HKMA20B/2008", "A/HongKong/HKMA20D/2008", "A/HongKong/HKMA20E/2008",
        "A/Busan/15453/2009", "A/Cambodia/NHRCC00001/2009", "A/Cambodia/NHRCC00002/2009",
        "A/Cambodia/NHRCC00003/2009", "A/Cambodia/NHRCC00006/2009", "A/Cambodia/NHRCC00007/2009",
        "A/Cambodia/NHRCC00008/2009", "A/Cambodia/NHRCC00010/2009", "A/Cambodia/NHRCC00011/2009",
        "A/Kansas/13/2009", "A/StPetersburg/5/2009", "A/Nanjing/1/2010", "A/Pennsylvania/14/2010",
        "A/Pennsylvania/40/2010", "A/Guam/AF2771/2011", "A/Indiana/8/2011", "A/Kenya/155/2011",
        "A/Kenya/168/2011", "A/Kenya/170/2011", "A/Nepal/142/2011", "A/Pennsylvania/09/2011",
        "A/Pennsylvania/9/2011", "A/Quebec/167936/2011", "A/Quebec/170658/2011", "A/India/6352/2012",
        "A/Indiana/7/2012", "A/Indiana/08/2012", "A/Indiana/8/2012", "A/Indiana/13/2012", "A/Ohio/18/2012",
        "A/Mendoza/176711/2012", "A/Ohio/34/2012", "A/Helsinki/942/2013", "A/Indiana/4/2013",
        "A/Indiana/5/2013", "A/Indiana/06/2013", "A/Indiana/6/2013", "A/Indiana/8/2013",
        "A/Indiana/11/2013", "A/Indiana/17/2013", "A/Indiana/21/2013", "A/Iowa/4/2013",
        "A/Jiangsu-Tianning/1707/2013", "A/HuNan/01/2014", "A/HuNan/1/2014", "A/Iran/140/2014",
        "A/Iran/227/2014", "A/Iran/234/2014", "A/Iran/407/2014", "A/Jiangsu-Chongchuan/1830/2014",
        "A/Jiangsu-Chongchuan/12179/2014", "A/Ohio/2/2014", "A/Ohio/4319/2014", "A/SaoPaulo/3-34891/2014",
        "A/Wisconsin/24/2014", "A/Jiangsu-Chongchuan/12714/2015", "A/Jiangsu-Chongchuan/12908/2015",
        "A/Michigan/39/2015", "A/NewJersey/53/2015", "A/SaoPaulo/36178/2015", "A/SaoPaulo/61282/2015",
        "A/SaoPaulo/65194/2015", "A/Sydney/53/2015", "A/Michigan/82/2016", "A/Michigan/83/2016",
        "A/Michigan/84/2016", "A/Michigan/87/2016", "A/Michigan/89/2016", "A/Michigan/90/2016",
        "A/Michigan/91/2016", "A/Michigan/93/2016", "A/Michigan/94/2016", "A/Michigan/95/2016",
        "A/Michigan/96/2016", "A/Ohio/27/2016", "A/Ohio/28/2016", "A/Ohio/32/2016", "A/Ohio/33/2016",
        "A/SouthAustralia/47/2016", "A/Ohio/35/2016", "A/Zhejiang-Wuxin/1300/2016",
        "A/London/16U363140-91_S43_L001/2016", "A/Catalonia/NSVH100533399/2017",
        "A/Catalonia/NSVH100560486/2017", "A/Piaui/494713/2017", "A/SaoPaulo/554097/2017",
        "A/Ireland/61451/2017", "A/Ireland/53051/2017", "A/Ireland/62518/2017", "A/Ireland/52590/2017",
        "A/Ireland/61097/2017", "A/Shanghai-Minxing/1482/2017"
    ],
    'h1n1pdm': [
        "A/Wisconsin/87/2005", "A/Iowa/1/2006", "A/Iowa/1/2006-egg", "A/Ohio/1/2007", "A/Ohio/2/2007",
        "A/Illinois/9/2007", "A/Minnesota/3/2008", "A/Texas/14/2008", "A/Malaysia/2142295/2009",
        "A/Malaysia/2142299/2009", "A/Singapore/SM15/2009", "A/HongKong/H090-665-V1/2009",
        "A/HongKong/H090-667-V1/2009", "A/HongKong/H090-751-V3/2009", "A/Austria/183/2009-egg",
        "A/HongKong/H090-771-V1/2009", "A/HongKong/H090-774-V1/2009", "A/HongKong/H090-684-V10/2009",
        "A/Austria/183/2009-egg", "A/Malaysia/2143696/2009", "A/HongKong/H090-667-V2/2009",
        "A/Shandong/1/2009-egg", "A/California/07/2009NYMC-X18113/198", "A/Iowa/2/2009",
        "A/Shandong/1/2009", "A/Iowa/2/2009-egg", "A/Asturias/RR6898/2010",
        "A/Christchurch/16/2010NIB-74xp13/202", "A/Tomsk/273-MA1/2010", "A/Wisconsin/28/2011",
        "A/Kenya/264/2012", "A/Missouri/12/2012", "A/Ontario/N163578/2012",
        "A/RioGrandedoNorte/117490/2012", "A/SriLanka/11/2012", "A/Arkansas/14/2013", "A/Arkansas/15/2013",
        "A/Hungary/02/2013", "A/Hungary/16/2013", "A/Poland/16/2013", "A/Brest/1161/2014", "A/Dakar/3/2014",
        "A/India/Pun1418633/2014", "A/Minnesota/33/2014", "A/Norway/1929/2014", "A/Iowa/39/2015",
        "A/India/4101/2015", "A/Minnesota/46/2015", "A/Ohio/9/2015", "A/Bari/166/2016", "A/Bari/167/2016",
        "A/Cherkessk/2/2016", "A/Belgium/G0027/2016"
    ],
    'vic':[
        "A/Malaysia/438/2016", "B/Togo/LNG/419/2013", "B/Brisbine/33/2008", "B/Kol/2024/2008",
        "B/Kolkata/1373/2008", "B/Kolkata/2024/2008", "B/Kolkata/372/2010", "B/Cambodia/26/2011",
        "B/Cambodia/30/2011", "B/Cambodia/62/2011", "B/Cambodia/89/2011", "B/Cambodia/V1005378/2011",
        "B/Stockholm/7/2011", "B/Bangkok/SI17/2012", "B/Bangkok/SI58/2012", "B/SouthAustralia/81/2012",
        "B/Netherlands/76/2014", "B/NewCaledonia/119/2015", "B/Thailand/CU-B11637/2015",
        "B/Brisbane/14/2016", "B/Sydney/6/2016", "B/Netherlands/883/2016"
    ],
    "yam":[
        "B/Kisumu/7/2005", "B/Nairobi/351/2005", "B/Kolkata/2546/2009", "B/Kolkata/N-1272/2009",
        "B/Kolkata/N-2047/2009", "B/Riyadh/3/2010", "B/Riyadh/4/2010", "B/England/581/2012",
        "B/Thailand/CU-B10303/2014", "B/Catalonia/NSVH100562319/2017", "B/Norway/2155/2017"
    ]
}

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
            "path": "metadata/h1n1pdm_outgroup.gb",
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
            "include": 0,
            "genes": ["SigPep", "HA1", "HA2"]
        },
        "na": {
            "path": "metadata/h1n1pdm_na_outgroup.gb",
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
        'A/Antananarivo/1067/2016', 'A/Florida/23/2017', 'A/Washington/16/2017', 'A/NorthCarolina/4/2017',
        'A/SouthCarolina/4/2017', 'A/Texas/71/2017', 'A/Wisconsin/327/2017', 'A/Texas/68/2017',
        'A/Brisbane/29/2017',
        'A/Texas/50/2012-egg', 'A/Switzerland/9715293/2013-egg', 'A/HongKong/7127/2014-egg',
        'A/NewCaledonia/71/2014-egg', 'A/HongKong/4801/2014-egg', 'A/Saitama/103/2014-egg',
        'A/Singapore/Infimh-16-0019/2016-egg', 'A/Norway/3806/2016-egg',
    ],
    'h1n1pdm':[
        'A/California/7/2009', 'A/Bangladesh/2021/2012', 'A/Victoria/367/2012', 'A/Brisbane/96/2012',
        'A/SouthAfrica/3626/2013', 'A/Brisbane/28/2013', 'A/Bolivia/559/2013', 'A/Florida/62/2014',
        'A/Newcaledonia/58/2014', 'A/Tasmania/24/2014', 'A/Florida/62/2014', 'A/Michigan/45/2015',
        'A/Perth/103/2015', 'A/Brisbane/181/2015', 'A/Singapore/GP1908/2015', 'A/Iowa/53/2015',
        'A/StPetersburg/61/2015', 'A/Minnesota/32/2015', 'A/Bangladesh/3002/2015', 'A/Indiana/21/2016',
        'A/Victoria/503/2016'
    ],
    'vic':[
        'B/Shangdong/7/1997', 'B/HongKong/330/2001', 'B/Malaysia/2506/2004', 'B/Ohio/1/2005',
        'B/Brisbane/60/2008', 'B/Utah/8/2012', 'B/Montana/5/2012', 'B/Texas/2/2013', 'B/Florida/33/2014',
        'B/Indiana/25/2015', 'B/Florida/78/2015', 'B/Colorado/6/2017'
    ],
    'yam':[
        'B/Beijing/184/1993', 'B/Sichuan/379/1999', 'B/Shanghai/361/2002', 'B/Florida/4/2006',
        'B/Wisconsin/1/2010', 'B/Massachusetts/2/2012', 'B/Phuket/3073/2013', 'B/Utah/9/2014',
        'B/Brisbane/9/2014', 'B/Guangdong-Liwan/1133/2014', 'B/California/12/2015', 'B/Arizona/10/2015',
        'B/NewHampshire/1/2016'
    ]
}

#Add list of (egg-passaged and matched) sequences to forcibly included
reference_viruses['h3n2'] = reference_viruses['h3n2'] + ['A/SaoPaulo/42047/2015-egg', 'A/Berlin/93/2011-egg', 'A/MatoGrossodoSul/559890/2017-egg', 'A/Corsica/15-05/2015-egg', 'A/SaoPaulo/10769/2014-egg', 'A/SaoPaulo/52727/2015-egg', 'A/Ohio/31/2017', 'A/Ohio/31/2017-egg', 'A/HongKong/4544/2009-egg', 'A/Singapore/GP2646/2016', 'A/Singapore/GP2646/2016-egg', 'A/Parana/490/2017-cell', 'A/Parana/490/2017-egg', 'A/Singapore/KK1100/2016', 'A/Singapore/KK1100/2016-egg', 'A/Colombo/159/2003-egg', 'A/Rondonia/18925/2015-egg', 'A/SaoPaulo/538363/2017-egg', 'A/Jiangsu-Donghai/57/2012-egg', 'A/Jiangsu-Chongchuan/12908/2015-egg', 'A/NewYork/2e/2006-egg', 'A/Perth/239/2013-egg', 'A/Kanagawa/AC1709/2017-cell', 'A/Kanagawa/AC1709/2017-egg', 'A/Romania/208407/2017-egg', 'A/SaoPaulo/452128/2017-egg', 'A/Malaysia/1199/2007-egg', 'A/Victoria/208/2009', 'A/Victoria/208/2009-cell', 'A/Victoria/208/2009-egg', 'A/SaoPaulo/6-10465/2014-egg', 'A/Washington/1e/2004-egg', 'A/Ontario/K4050772/2010-egg', 'A/SouthAfrica/9015/2013-egg', 'A/Antananarivo/1067/2016', 'A/Antananarivo/1067/2016-cell', 'A/Antananarivo/1067/2016-egg', 'A/HongKong/1-5/1968-egg', 'A/Ontario/C720531/2010-egg', 'A/Parana/717/2017-egg', 'A/MatoGrosso/54934/2015-egg', 'A/HoChiMinh/6077/2010-egg', 'A/CzechRepublic/28/2015-egg', 'A/SaoPaulo/594145/2017-egg', 'A/HongKong/26560/2009', 'A/HongKong/26560/2009-cell', 'A/HongKong/26560/2009-egg', 'A/Illinois/14/2008-cell', 'A/Illinois/14/2008-egg', 'A/Indiana/38/2012-cell', 'A/Indiana/38/2012-egg', 'A/Singapore/H2011479/2011-cell', 'A/Singapore/H2011479/2011-egg', 'A/Singapore/C2011173/2011-cell', 'A/Singapore/C2011173/2011-egg', 'A/Beijing/32/1992', 'A/Beijing/32/1992-cell', 'A/Beijing/32/1992-egg', 'A/HongKong/1-1-MA-12B/1968-egg', 'A/Taiwan/1021/1982-egg', 'A/SaoPaulo/6-14388/2015-egg', 'A/Newcastle/25/2014-egg', 'A/SaoPaulo/53072/2014-egg', 'A/SaoPaulo/340809/2016-egg', 'A/Brazil/2363/2005-egg', 'A/Switzerland/7586983/2012-cell', 'A/Switzerland/7586983/2012-egg', 'A/Hawaii/16e/2004-egg', 'A/Sydney/103/2017-cell', 'A/Sydney/103/2017-egg', 'A/SouthAfrica/6451/2013-egg', 'A/Brisbane/47/2015-cell', 'A/Brisbane/47/2015-egg', 'A/CostaRica/4700/2013-cell', 'A/CostaRica/4700/2013-egg', 'A/SaoPaulo/2-16687/2015-egg', 'A/SaoPaulo/4-12847/2014-egg', 'A/Quebec/26-310106/2006-egg', 'A/Ontario/C687409/2010-egg', 'A/CzechRepublic/33/2015-cell', 'A/CzechRepublic/33/2015-egg', 'A/Corsica/3521-06/2015-egg', 'A/Singapore/H2011704C/2011-egg', 'A/Glasgow/RVL5/2012-egg', 'A/SaoPaulo/51790/2015-egg', 'A/Ontario/C575478/2010-egg', 'A/Newcaledonia/71/2014-cell', 'A/Newcaledonia/71/2014-egg', 'A/SaoPaulo/549896/2017-egg', 'A/SaoPaulo/535746/2017-egg', 'A/Fukuoka/55/2002', 'A/Fukuoka/55/2002-egg', 'A/UnitedKingdom/1864/2003-egg', 'A/Bilthoven/808/1969', 'A/Bilthoven/808/1969-egg', 'A/HongKong/4593e/2004-egg', 'A/Beijing/353/1989', 'A/Beijing/353/1989-cell', 'A/Beijing/353/1989-egg', 'A/SaoPaulo/3-24388/2014-egg', 'A/Beijing/352/1989', 'A/Beijing/352/1989-egg', 'A/Ontario/C582144/2010-egg', 'A/Kentucky/6e/2004-egg', 'A/Ohio/2/2012-cell', 'A/Ohio/2/2012-egg', 'A/Parana/484/2017-egg', 'A/Philippines/472/2002', 'A/Philippines/472/2002-cell', 'A/Philippines/472/2002-egg', 'A/SouthDakota/16/2012-egg', 'A/Parana/722/2017-egg', 'A/Goias/359849/2017-egg', 'A/Kanagawa/IC1739/2018-cell', 'A/Kanagawa/IC1739/2018-egg', 'A/Bangkok/1/1997', 'A/Bangkok/1/1997-egg', 'A/Singapore/H2013721d/2013-egg', 'A/Fiji/60/2016', 'A/Fiji/60/2016-cell', 'A/Fiji/60/2016-egg', 'A/Hawaii/23/2012-cell', 'A/Hawaii/23/2012-egg', 'A/SaoPaulo/11151/2014-egg', 'A/Ontario/C627683/2010-egg', 'A/HongKong/1-1-MA-20D/1968-egg', 'A/HoChiMinh/5517/2010-egg', 'A/SaoPaulo/510403/2017-egg', 'A/DistritoFederal/58903/2015-egg', 'A/Quebec/26-041105/2005-cell', 'A/Quebec/26-041105/2005-egg', 'A/Romania/3781-3820/2016-egg', 'A/HongKong/1-1-MA-12C/1968-egg', 'A/MinasGerais/625/2017-egg', 'A/HongKong/7127/2014-cell', 'A/HongKong/7127/2014-egg', 'A/Singapore/H2010211/2010-cell', 'A/Singapore/H2010211/2010-egg', 'A/Singapore/H2011570/2011-cell', 'A/Singapore/H2011570/2011-egg', 'A/Singapore/C2011458/2011-cell', 'A/Singapore/C2011458/2011-egg', 'A/Hanoi/Q148/2007-egg', 'A/Maine/5/2012-cell', 'A/Maine/5/2012-egg', 'A/Ontario/C706265/2010-egg', 'A/Kentucky/5/2011-cell', 'A/Kentucky/5/2011-egg', 'A/Catalonia/585566S/2016-egg', 'A/HongKong/1052/2003', 'A/HongKong/1052/2003-cell', 'A/HongKong/1052/2003-egg', 'A/Panama/2007/1999', 'A/Panama/2007/1999-cell', 'A/Panama/2007/1999-egg', 'A/Singapore/H2011482/2011-cell', 'A/Singapore/H2011482/2011-egg', 'A/SaoPaulo/6-14390/2015-egg', 'A/Glasgow/407664/2012', 'A/Glasgow/407664/2012-cell', 'A/Glasgow/407664/2012-egg', 'A/Goias/53024/2014-egg', 'A/Singapore/H2011471c/2011-cell', 'A/Singapore/H2011471c/2011-egg', 'A/Victoria/505/2013-cell', 'A/Victoria/505/2013-egg', 'A/Rondonia/50872/2013-egg', 'A/HongKong/1-5-MA21-1/1968', 'A/HongKong/1-5-MA21-1/1968-egg', 'A/SaoPaulo/3-37389/2015-egg', 'A/Perth/4/2008', 'A/Perth/4/2008-egg', 'A/CzechRepublic/76/2012-egg', 'A/Singapore/H2009485/2009-cell', 'A/Singapore/H2009485/2009-egg', 'A/SaoPaulo/9-134/2014-egg', 'A/Ontario/C1995/2010-egg', 'A/SaoPaulo/1054970/2017-egg', 'A/Tocantins/24996/2015-egg', 'A/Norway/4293/2016-cell', 'A/Norway/4293/2016-egg', 'A/NewYork/3e/2006-cell', 'A/NewYork/3e/2006-egg', 'A/Victoria/5006/2014-cell', 'A/Victoria/5006/2014-egg', 'A/Albany/4/1969-egg', 'A/Singapore/GP2050/2015-cell', 'A/Singapore/GP2050/2015-egg', 'A/MatoGrosso/350652/2016-egg', 'A/Albany/18/1968-egg', 'A/Honduras/3112/2006-cell', 'A/Honduras/3112/2006-egg', 'A/Singapore/C2011512/2011-cell', 'A/Singapore/C2011512/2011-egg', 'A/Singapore/C2011825/2011-cell', 'A/Singapore/C2011825/2011-egg', 'A/Indiana/60/2012-cell', 'A/Indiana/60/2012-egg', 'A/Kanagawa/IC1619/2017-cell', 'A/Kanagawa/IC1619/2017-egg', 'A/HongKong/1-4-MA21-1/1968', 'A/HongKong/1-4-MA21-1/1968-egg', 'A/SaoPaulo/3-32002/2015-egg', 'A/Missouri/4/2013-cell', 'A/Missouri/4/2013-egg', 'A/California/32/1999-egg', 'A/NewYork/39/2012-cell', 'A/NewYork/39/2012-egg', 'A/Arizona/62e/2003-egg', 'A/SriLanka/61/2015-cell', 'A/SriLanka/61/2015-egg', 'A/Massachusetts/11/2013-cell', 'A/Massachusetts/11/2013-egg', 'A/Memphis/11/1986-egg', 'A/Corsica/15-12/2015-egg', 'A/Victoria/361/2011-cell', 'A/Victoria/361/2011-egg', 'A/Victoria/6001/2015-egg', 'A/SaoPaulo/34136/2014-egg', 'A/Albany/17/1968-egg', 'A/Jiangsu-Wujin/318/2013-egg', 'A/Jilin/18/2003', 'A/Jilin/18/2003-egg', 'A/Johannesburg/31/2000-egg', 'A/GuangdongFutian/1253/2009', 'A/GuangdongFutian/1253/2009-egg', 'A/Ontario/P21548/2010-egg', 'A/Wisconsin/24/2014-cell', 'A/Wisconsin/24/2014-egg', 'A/Kanagawa/ZC1710/2017-cell', 'A/Kanagawa/ZC1710/2017-egg', 'A/Victoria/605/1998-egg', 'A/Parana/725/2017-egg', 'A/Hawaii/15e/2004-egg', 'A/Tocantins/37342/2015-egg', 'A/Kentucky/24/2012-egg', 'A/Corsica/3975-10/2015-egg', 'A/Singapore/H2013060/2013-egg', 'A/Wisconsin/67e5/2005-egg', 'A/SaoPaulo/20768/2015-egg', 'A/Goias/6394/2015-egg', 'A/Mexico/InDRE2112/2005-egg', 'A/Singapore/C2011518/2011-cell', 'A/Singapore/C2011518/2011-egg', 'A/Bilthoven/10684/1982', 'A/Bilthoven/10684/1982-egg', 'A/Victoria/710/2013', 'A/Victoria/710/2013-egg', 'A/RhodeIsland/1/2018', 'A/RhodeIsland/1/2018-egg', 'A/Parana/489/2017-egg', 'A/Stockholm/36/2012', 'A/Stockholm/36/2012-cell', 'A/Stockholm/36/2012-egg', 'A/Florida/40/2012-egg', 'A/Memphis/1/1971', 'A/Memphis/1/1971-egg', 'A/Singapore/H2013901/2013-egg', 'A/HoChiMinh/6328/2010-egg', 'A/SaoPaulo/3-24386/2014-egg', 'A/Fiji/1/2014-cell', 'A/Fiji/1/2014-egg', 'A/Jiangsu-Danyang/1647/2012-egg', 'A/SaoPaulo/593870/2017-egg', 'A/Japan/120/2004', 'A/Japan/120/2004-egg', 'A/SaoPaulo/9-1220/2015-egg', 'A/Perth/201/2001-cell', 'A/Perth/201/2001-egg', 'A/SouthAfrica/5853/2013-egg', 'A/Singapore/C2009515/2009-cell', 'A/Singapore/C2009515/2009-egg', 'A/Norway/2605/2015-cell', 'A/Norway/2605/2015-egg', 'A/SaoPaulo/706496/2018-egg', 'A/CzechRepublic/40/2015-egg', 'A/Wellington/1/2004', 'A/Wellington/1/2004-cell', 'A/Wellington/1/2004-egg', 'A/HongKong/2142/2009-cell', 'A/HongKong/2142/2009-egg', 'A/California/61/2013-cell', 'A/California/61/2013-egg', 'A/SouthAfrica/9000/2013-egg', 'A/Micronesia/6689/2012-egg', 'A/Singapore/37/2009-egg', 'A/Singapore/37/2004', 'A/Singapore/37/2004-egg', 'A/CzechRepublic/700/2014-egg', 'A/Udorn/1972', 'A/Udorn/1972-egg', 'A/Taiwan/613/2012-egg', 'A/CzechRepublic/39/2015-egg', 'A/SouthDakota/15/2012-egg', 'A/SaoPaulo/614558/2017-egg', 'A/SaoPaulo/557516/2017-egg', 'A/Egypt/130/2002-egg', 'A/DistritoFederal/89664/2014-egg', 'A/Brisbane/1/2005-egg', 'A/Singapore/H2009518/2009-cell', 'A/Singapore/H2009518/2009-egg', 'A/Brisbane/4/2017-cell', 'A/Brisbane/4/2017-egg', 'A/Philippines/2/1982-cell', 'A/Philippines/2/1982-egg', 'A/SaoPaulo/2933/2014-egg', 'A/Kentucky/3/2018', 'A/Kentucky/3/2018-egg', 'A/CzechRepublic/119/2012-egg', 'A/CostaRica/4009/2013-cell', 'A/CostaRica/4009/2013-egg', 'A/Romania/209119/2017-egg', 'A/Montana/3e/2005-egg', 'A/Hawaii/22/2012-cell', 'A/Hawaii/22/2012-egg', 'A/SaoPaulo/47588/2015-egg', 'A/Alaska/23/2012-cell', 'A/Alaska/23/2012-egg', 'A/Sydney/195/2012-cell', 'A/Sydney/195/2012-egg', 'A/UlanUde/1/2000-egg', 'A/Goias/386531/2017-egg', 'A/Victoria/513/2004', 'A/Victoria/513/2004-cell', 'A/Victoria/513/2004-egg', 'A/SaoPaulo/54274/2014-egg', 'A/CzechRepublic/23/2014-egg', 'A/SaoPaulo/24295/2015-egg', 'A/Brisbane/3/2005', 'A/Brisbane/3/2005-cell', 'A/Brisbane/3/2005-egg', 'A/SouthAfrica/7241/2013-egg', 'A/Corsica/40-04/2015-egg', 'A/Lisbon/SU91/2012-cell', 'A/Lisbon/SU91/2012-egg', 'A/SaoPaulo/5-13834/2014-egg', 'A/Canberra/2/2014-cell', 'A/Canberra/2/2014-egg', 'A/Singapore/H2011515/2011-cell', 'A/Singapore/H2011515/2011-egg', 'A/SantaCatarina/452/2017-egg', 'A/Jilin/14/2003', 'A/Jilin/14/2003-egg', 'A/SaoPaulo/510173/2017-egg', 'A/HongKong/2277/2017', 'A/HongKong/2277/2017-cell', 'A/HongKong/2277/2017-egg', 'A/Georgia/169/2017-egg', 'A/SaoPaulo/13-8487/2015-egg', 'A/CzechRepublic/44/2015-egg', 'A/Memphis/102/1972', 'A/Memphis/102/1972-egg', 'A/SaoPaulo/51373/2014-egg', 'A/Indiana/5/2013', 'A/Indiana/5/2013-cell', 'A/Indiana/5/2013-egg', 'A/MatoGrossodoSul/62052/2015-egg', 'A/SaoPaulo/3-29553/2014-egg', 'A/Singapore/CGH05/2013-cell', 'A/Singapore/CGH05/2013-egg', 'A/Canberra/3/2014-egg', 'A/Missouri/13e/2005-egg', 'A/Missouri/13e/2004-egg', 'A/HongKong/46/1980-egg', 'A/Sichuan/2/1987', 'A/Sichuan/2/1987-egg', 'A/Wisconsin/19/2017', 'A/Wisconsin/19/2017-cell', 'A/Wisconsin/19/2017-egg', 'A/SaoPaulo/13-8806/2015-egg', 'A/Singapore/GP1523/2014-egg', 'A/SaoPaulo/9-10272/2015-egg', 'A/Texas/40e/2003-egg', 'A/Goias/375493/2017-egg', 'A/SaoPaulo/6-11670/2014-egg', 'A/Utah/10/2012', 'A/Utah/10/2012-egg', 'A/Quebec/26-031205/2005-cell', 'A/Quebec/26-031205/2005-egg', 'A/SaoPaulo/7-5547/2014-egg', 'A/Romania/208404/2017-egg', 'A/Stockholm/6/2014', 'A/Stockholm/6/2014-cell', 'A/Stockholm/6/2014-egg', 'A/Christchurch/28/2003', 'A/Christchurch/28/2003-cell', 'A/Christchurch/28/2003-egg', 'A/MatoGrossodoSul/70736/2015-egg', 'A/SaoPaulo/56468/2014-egg', 'A/Brisbane/210/2010', 'A/Brisbane/210/2010-cell', 'A/Brisbane/210/2010-egg', 'A/HongKong/1-10-MA21-3/1968', 'A/HongKong/1-10-MA21-3/1968-egg', 'A/Brisbane/340/2003', 'A/Brisbane/340/2003-egg', 'A/SaoPaulo/48012/2015-egg', 'A/CzechRepublic/126/2012-egg', 'A/Ontario/K4023106/2010-egg', 'A/SouthDakota/12/2012-cell', 'A/SouthDakota/12/2012-egg', 'A/EspiritoSanto/500/2017-egg', 'A/Perth/10/2010', 'A/Perth/10/2010-cell', 'A/Perth/10/2010-egg', 'A/Wisconsin/15/2009', 'A/Wisconsin/15/2009-cell', 'A/Wisconsin/15/2009-egg', 'A/SantaCatarina/454/2017-egg', 'A/Singapore/H2010471C/2010-egg', 'A/SaoPaulo/662138/2017-egg', 'A/Mississippi/4/2008-egg', 'A/Brisbane/24/2008', 'A/Brisbane/24/2008-cell', 'A/Brisbane/24/2008-egg', 'A/SaoPaulo/4-21208/2015-egg', 'A/Washington/42/2007-egg', 'A/Tunisia/723/2009-egg', 'A/Quebec/26-151105/2005-egg', 'A/SaoPaulo/11-29688/2013-egg', 'A/CzechRepublic/30/2015-egg', 'A/Ontario/C720525/2010-egg', 'A/Singapore/H2010822/2010-cell', 'A/Singapore/H2010822/2010-egg', 'A/Tasmania/18/2017-cell', 'A/Tasmania/18/2017-egg', 'A/Singapore/H2011474/2011-cell', 'A/Singapore/H2011474/2011-egg', 'A/SaoPaulo/6-11921/2014-egg', 'A/RhodeIsland/2/2018', 'A/RhodeIsland/2/2018-egg', 'A/Singapore/H2013403/2013-egg', 'A/CzechRepublic/131/2012-egg', 'A/Kanagawa/AC1704/2017-egg', 'A/SaoPaulo/733190/2018-egg', 'A/Kanagawa/IC1610/2017-cell', 'A/Kanagawa/IC1610/2017-egg', 'A/RioGrandedoSul/661/2017-egg', 'A/DistritoFederal/58900/2015-egg', 'A/ChongqingYuzhong/169/2007', 'A/ChongqingYuzhong/169/2007-egg', 'A/Singapore/H2010389/2010-cell', 'A/Singapore/H2010389/2010-egg', 'A/Ontario/U10611/2010-egg', 'A/SaoPaulo/52081/2015-egg', 'A/NewHampshire/11/2012-cell', 'A/NewHampshire/11/2012-egg', 'A/Singapore/H2011447/2011-cell', 'A/Singapore/H2011447/2011-egg', 'A/SaoPaulo/548298/2017-egg', 'A/SaoPaulo/549098/2017-egg', 'A/Wuhan/359/1995', 'A/Wuhan/359/1995-egg', 'A/Alaska/4/2012-cell', 'A/Alaska/4/2012-egg', 'A/Singapore/TT1374/2016', 'A/Singapore/TT1374/2016-egg', 'A/Ontario/C835944/2010-egg', 'A/Alaska/5/2010-cell', 'A/Alaska/5/2010-egg', 'A/Kentucky/1/2018', 'A/Kentucky/1/2018-egg', 'A/Idaho/33/2016', 'A/Idaho/33/2016-cell', 'A/Idaho/33/2016-egg', 'A/Mississippi/5e/2004-egg', 'A/California/20e/2005-egg', 'A/Singapore/C2011384/2011-cell', 'A/Singapore/C2011384/2011-egg', 'A/Corsica/3975-04/2014-egg', 'A/SaoPaulo/596791/2017-egg', 'A/SaoPaulo/10-23957/2015-egg', 'A/Singapore/H2011471a/2011-cell', 'A/Singapore/H2011471a/2011-egg', 'A/Florida/21/2011-cell', 'A/Florida/21/2011-egg', 'A/Brisbane/2/2008', 'A/Brisbane/2/2008-egg', 'A/MinasGerais/66/2017', 'A/MinasGerais/66/2017-egg', 'A/Maine/2/2018', 'A/Maine/2/2018-egg', 'A/Singapore/H2010666C/2010-egg', 'A/Palau/6759/2014-cell', 'A/Palau/6759/2014-egg', 'A/Texas/71/2017', 'A/Texas/71/2017-cell', 'A/Texas/71/2017-egg', 'A/Sendai-H/F131/2006-egg', 'A/Corsica/40-03/2015-egg', 'A/Wisconsin/21e/2004-egg', 'A/Singapore/H2009679/2009-cell', 'A/Singapore/H2009679/2009-egg', 'A/Wisconsin/3/2007-cell', 'A/Wisconsin/3/2007-egg', 'A/Goias/49684/2015-egg', 'A/Switzerland/9715293/2013-cell', 'A/Switzerland/9715293/2013-egg', 'A/Stockholm/40/2012', 'A/Stockholm/40/2012-cell', 'A/Stockholm/40/2012-egg', 'A/Ontario/C37332/2011-egg', 'A/Goias/3910/2014-egg', 'A/SaoPaulo/6-6304/2015-egg', 'A/SaoPaulo/3-22794/2015-egg', 'A/Peru/5107/2006', 'A/Peru/5107/2006-egg', 'A/SouthAfrica/D1147/1996-egg', 'A/Wisconsin/16e/2004-egg', 'A/SaoPaulo/3-29379/2015-egg', 'A/Brisbane/342/2003', 'A/Brisbane/342/2003-cell', 'A/Brisbane/342/2003-egg', 'A/SaoPaulo/10-22878/2015-egg', 'A/Singapore/H2011499a/2011-cell', 'A/Singapore/H2011499a/2011-egg', 'A/MinasGerais/626/2017-egg', 'A/Singapore/H2010559C/2010-egg', 'A/Singapore/H2011463/2011-cell', 'A/Singapore/H2011463/2011-egg', 'A/Philippines/2191/2009', 'A/Philippines/2191/2009-cell', 'A/Philippines/2191/2009-egg', 'A/Palau/6775/2014-cell', 'A/Palau/6775/2014-egg', 'A/Singapore/C2011641/2011-cell', 'A/Singapore/C2011641/2011-egg', 'A/Victoria/500/2008', 'A/Victoria/500/2008-egg', 'A/Romania/206417/2017-egg', 'A/Ontario/C582093/2010-egg', 'A/CzechRepublic/34/2015-egg', 'A/SantaCatarina/479/2017-egg', 'A/Parana/491/2017-egg', 'A/SaoPaulo/10898/2014-egg', 'A/Corsica/4069-20/2015-egg', 'A/Oklahoma/5/2013-cell', 'A/Oklahoma/5/2013-egg', 'A/Taiwan/839/2009', 'A/Taiwan/839/2009-cell', 'A/Taiwan/839/2009-egg', 'A/HongKong/1-8-MA21-2/1968-egg', 'A/Peru/5102/2006', 'A/Peru/5102/2006-egg', 'A/SaoPaulo/3-24357/2014-egg', 'A/Parana/719/2017-egg', 'A/NorthernTerritory/60-JY2/1968-egg', 'A/SaoPaulo/2-1030/2015-egg', 'A/Romania/225758/2018-egg', 'A/Greece/4/2017-cell', 'A/Greece/4/2017-egg', 'A/Norway/4465/2016', 'A/Norway/4465/2016-cell', 'A/Norway/4465/2016-egg', 'A/Bangladesh/3012/2009', 'A/Bangladesh/3012/2009-egg', 'A/California/7/2004', 'A/California/7/2004-cell', 'A/California/7/2004-egg', 'A/SaoPaulo/539114/2017-egg', 'A/Alaska/3e/2004-egg', 'A/Brisbane/82/2015-cell', 'A/Brisbane/82/2015-egg', 'A/SouthAustralia/274/2017-cell', 'A/SouthAustralia/274/2017-egg', 'A/Johannesburg/92/1999-egg', 'A/Corsica/15-06/2015-egg', 'A/Colorado/18/2012-cell', 'A/Colorado/18/2012-egg', 'A/Singapore/INFTT-16-0612/2016', 'A/Singapore/INFTT-16-0612/2016-egg', 'A/HoChiMinh/5897/2010-egg', 'A/Brisbane/318/2016-cell', 'A/Brisbane/318/2016-egg', 'A/Denmark/48/2014-cell', 'A/Denmark/48/2014-egg', 'A/SaoPaulo/10-20095/2015-egg', 'A/Victoria/512/2004-egg', 'A/SaoPaulo/3308/2014-egg', 'A/Kanagawa/IC1725/2018-cell', 'A/Kanagawa/IC1725/2018-egg', 'A/Memphis/109/1972-egg', 'A/Georgia/238/2017-egg', 'A/HoChiMinh/85410/2010-egg', 'A/NewYork/5/2014-cell', 'A/NewYork/5/2014-egg', 'A/CzechRepublic/94/2015-egg', 'A/Wisconsin/67/2005', 'A/Wisconsin/67/2005-cell', 'A/Wisconsin/67/2005-egg', 'A/Tasmania/11/2014-cell', 'A/Tasmania/11/2014-egg', 'A/HongKong/1-1-MA-12E/1968-egg', 'A/CzechRepublic/36/2015-egg', 'A/SouthAfrica/9162/2013-egg', 'A/Singapore/H2011460cC/2011-egg', 'A/Wyoming/25e/2005-egg', 'A/SouthAfrica/8195/2013-egg', 'A/HoChiMinh/6127/2010-egg', 'A/Victoria/506/2013-cell', 'A/Victoria/506/2013-egg', 'A/Moscow/10/1999', 'A/Moscow/10/1999-egg', 'A/Goias/44080/2015-egg', 'A/SaoPaulo/78870/2015-egg', 'A/SaoPaulo/6-726/2015-egg', 'A/CzechRepublic/21/2014-egg', 'A/Shanghai/11/1987', 'A/Shanghai/11/1987-cell', 'A/Shanghai/11/1987-egg', 'A/HoChiMinh/6528/2010-egg', 'A/SaoPaulo/38433/2015-egg', 'A/Romania/206324/2017-egg', 'A/Southaustralia/55/2014-cell', 'A/Southaustralia/55/2014-egg', 'A/Romania/206784/2017-egg', 'A/Wyoming/3e/2003-egg', 'A/Indiana/21/2012-cell', 'A/Indiana/21/2012-egg', 'A/Mexico/29/2006-egg', 'A/Brazil/2353/2005-egg', 'A/NewJersey/2e/2006-egg', 'A/Colombo/190/2003-egg', 'A/Memphis/7/1985-egg', 'A/Memphis/8/1988-egg', 'A/Ontario/C4144/2011-egg', 'A/Pennsylvania/CHOP10/2012-cell', 'A/Pennsylvania/CHOP10/2012-egg', 'A/DistritoFederal/51671/2015-egg', 'A/Ontario/N91738/2010-egg', 'A/Shandong-Laicheng/1763/2016-egg', 'A/Romania/207793/2017-egg', 'A/Texas/1/1977', 'A/Texas/1/1977-cell', 'A/Texas/1/1977-egg', 'A/Almaty/2958/2013-cell', 'A/Almaty/2958/2013-egg', 'A/WestVirginia/6/2011', 'A/WestVirginia/6/2011-cell', 'A/WestVirginia/6/2011-egg', 'A/Aichi/2/1968', 'A/Aichi/2/1968-egg', 'A/DistritoFederal/75551/2014-egg', 'A/SaoPaulo/348596/2016-egg', 'A/Taiwan/1363/2012-egg', 'A/SaoPaulo/36178/2015-egg', 'A/Singapore/22/2012-cell', 'A/Singapore/22/2012-egg', 'A/Parana/481/2017-egg', 'A/Romania/3748-3722/2016-egg', 'A/SouthAfrica/8210/2013-egg', 'A/Romania/225833/2018-cell', 'A/Romania/225833/2018-egg', 'A/HongKong/1-1-MA-12A/1968-egg', 'A/CzechRepublic/42/2015-egg', 'A/Singapore/H2011504C/2011-egg', 'A/Jiangsu-Xinpu/1205/2012-egg', 'A/HongKong/45/1980-egg', 'A/Alaska/232/2015', 'A/Alaska/232/2015-cell', 'A/Alaska/232/2015-egg', 'A/NorthCarolina/13/2014-cell', 'A/NorthCarolina/13/2014-egg', 'A/Pennsylvania/CHOP1/2013-cell', 'A/Pennsylvania/CHOP1/2013-egg', 'A/Missouri/2/2012-cell', 'A/Missouri/2/2012-egg', 'A/SaoPaulo/5-12085/2015-egg', 'A/Bangladesh/1817/2005-egg', 'A/Brisbane/183/2015-cell', 'A/Brisbane/183/2015-egg', 'A/SaoPaulo/3-36419/2014-egg', 'A/HongKong/2286/2017', 'A/HongKong/2286/2017-cell', 'A/HongKong/2286/2017-egg', 'A/SaoPaulo/3-35437/2015-egg', 'A/SaoPaulo/2932/2014-egg', 'A/Corsica/15-08/2015-egg', 'A/Alaska/240/2015', 'A/Alaska/240/2015-egg', 'A/Acores/SU43/2012-cell', 'A/Acores/SU43/2012-egg', 'A/Iowa/7/2011', 'A/Iowa/7/2011-cell', 'A/Iowa/7/2011-egg', 'A/CzechRepublic/29/2015-egg', 'A/Brisbane/204/2013-egg', 'A/SaoPaulo/610485/2017-egg', 'A/Arizona/13/2010-cell', 'A/Arizona/13/2010-egg', 'A/NewYork/55/2004', 'A/NewYork/55/2004-egg', 'A/HongKong/1-6-MA21-3/1968-egg', 'A/Romania/206972/2017-egg', 'A/Saopaulo/603559/2017-egg', 'A/Goias/387153/2017-egg', 'A/SaoPaulo/14332/2014-egg', 'A/Sydney/37/2014-cell', 'A/Sydney/37/2014-egg', 'A/Taiwan/285/2012-egg', 'A/MatoGrosso/45246/2014-egg', 'A/Singapore/C2009485a/2009-cell', 'A/Singapore/C2009485a/2009-egg', 'A/Romania/207464/2017-egg', 'A/Tunisia/1851/2007-egg', 'A/Auckland/10/2015-cell', 'A/Auckland/10/2015-egg', 'A/CzechRepublic/31/2015-egg', 'A/NorthernTerritory/60/1968', 'A/NorthernTerritory/60/1968-egg', 'A/Pennsylvania/CHOP6/2013-cell', 'A/Pennsylvania/CHOP6/2013-egg', 'A/Victoria/673/2014-cell', 'A/Victoria/673/2014-egg', 'A/Goias/358931/2017-egg', 'A/SaoPaulo/2-16355/2015-egg', 'A/SouthAustralia/24/2014-egg', 'A/EspiritoSanto/502/2017-egg', 'A/Ontario/C610464/2010-egg', 'A/SaoPaulo/911967/2017-egg', 'A/Piaui/32789/2015-egg', 'A/Guadeloupe/202/2010', 'A/Guadeloupe/202/2010-cell', 'A/Guadeloupe/202/2010-egg', 'A/SouthAfrica/7073/2013-egg', 'A/Jiangsu-Chongchuan/1830/2014-egg', 'A/SaoPaulo/50168/2015-egg', 'A/SouthAfrica/2982/2015-cell', 'A/SouthAfrica/2982/2015-egg', 'A/Goias/388164/2017-egg', 'A/Netherlands/363/2006-egg', 'A/Singapore/C2011493/2011-cell', 'A/Singapore/C2011493/2011-egg', 'A/Shandong/9/1993', 'A/Shandong/9/1993-cell', 'A/Shandong/9/1993-egg', 'A/SaoPaulo/486371/2017-egg', 'A/SaoPaulo/6-11780/2014-egg', 'A/CzechRepublic/138/2012-egg', 'A/HongKong/4494/2009-egg', 'A/AthensGR/14/2013', 'A/AthensGR/14/2013-cell', 'A/AthensGR/14/2013-egg', 'A/Corsica/3975-03/2014-egg', 'A/Ontario/C706264/2010-egg', 'A/NewYork/40e/2004-egg', 'A/SaoPaulo/2-18547/2015-egg', 'A/Piaui/507920/2017-egg', 'A/Memphis/13/1988-egg', 'A/HongKong/5738/2014-cell', 'A/HongKong/5738/2014-egg', 'A/SaoPaulo/10770/2014-egg', 'A/Nevada/9/2014-cell', 'A/Nevada/9/2014-egg', 'A/SaoPaulo/10-23958/2015-egg', 'A/Geneva/1479911/2009', 'A/Geneva/1479911/2009-egg', 'A/Pennsylvania/11/2012-cell', 'A/Pennsylvania/11/2012-egg', 'A/Corsica/4076-04/2015-egg', 'A/Kanagawa/AC1620/2017-cell', 'A/Kanagawa/AC1620/2017-egg', 'A/Brisbane/71/2015-cell', 'A/Brisbane/71/2015-egg', 'A/HongKong/1-4-MA21-3/1968-egg', 'A/Perth/1/2014-egg', 'A/Wisconsin/19e/2004-egg', 'A/HoChiMinh/6268/2010-egg', 'A/SaoPaulo/34233/2014-egg', 'A/NewYork/38/2012-cell', 'A/NewYork/38/2012-egg', 'A/SantaCatarina/449/2017-egg', 'A/Goias/42119/2014-egg', 'A/SaoPaulo/3-36025/2015-egg', 'A/SaoPaulo/13-10842/2014-egg', 'A/SaoPaulo/12-10904/2015-egg', 'A/SaoPaulo/53924/2014-egg', 'A/HongKong/1-4/1968-egg', 'A/Singapore/C2011471/2011-cell', 'A/Singapore/C2011471/2011-egg', 'A/Singapore/Inftt-16-0612/2016-egg', 'A/SaoPaulo/54116/2014-egg', 'A/Singapore/C2011564/2011-cell', 'A/Singapore/C2011564/2011-egg', 'A/Taiwan/2361/2007-egg', 'A/NorthDakota/1e/2004-egg', 'A/Singapore/GP1599/2013-egg', 'A/Rondonia/53925/2015-egg', 'A/MatoGrossodoSul/523791/2017-egg', 'A/Kanagawa/IC1624/2017', 'A/Kanagawa/IC1624/2017-cell', 'A/Kanagawa/IC1624/2017-egg', 'A/Singapore/H2010797/2010-cell', 'A/Singapore/H2010797/2010-egg', 'A/Brisbane/9/2006', 'A/Brisbane/9/2006-cell', 'A/Brisbane/9/2006-egg', 'A/Harbin/15/1992', 'A/Harbin/15/1992-egg', 'A/Idaho/37/2016', 'A/Idaho/37/2016-cell', 'A/Idaho/37/2016-egg', 'A/Singapore/H2009510/2009-cell', 'A/Singapore/H2009510/2009-egg', 'A/Missouri/5/2008-cell', 'A/Missouri/5/2008-egg', 'A/Ontario/C835943/2010-egg', 'A/SaoPaulo/12-15655/2014-egg', 'A/Washington/13/2013-cell', 'A/Washington/13/2013-egg', 'A/HoChiMinh/6608/2010-egg', 'A/Ontario/C698767/2010-egg', 'A/Victoria/511/2015-cell', 'A/Victoria/511/2015-egg', 'A/SouthAustralia/40/2014-egg', 'A/Indiana/13/2012', 'A/Indiana/13/2012-egg', 'A/Norway/4849/2016-cell', 'A/Norway/4849/2016-egg', 'A/Ohio/13/2012', 'A/Ohio/13/2012-egg', 'A/SaoPaulo/88083/2013-egg', 'A/Mexico/InDRE29/2006-egg', 'A/HoChiMinh/7228/2010-egg', 'A/Kansas/14/2017', 'A/Kansas/14/2017-cell', 'A/Kansas/14/2017-egg', 'A/Goias/32286/2015-egg', 'A/Pennsylvania/CHOP9/2013-cell', 'A/Pennsylvania/CHOP9/2013-egg', 'A/Victoria/JY2/1968-egg', 'A/ShanghaiLuwan/188/2011', 'A/ShanghaiLuwan/188/2011-egg', 'A/HongKong/50/2016-cell', 'A/HongKong/50/2016-egg', 'A/SouthAfrica/8208/2013-egg', 'A/Romania/206613/2017-egg', 'A/Romania/206783/2017-egg', 'A/Hubei-Wuchang/170/2014-egg', 'A/Fiji/7/2015-cell', 'A/Fiji/7/2015-egg', 'A/Singapore/H2010564C/2010-egg', 'A/Ohio/1/2013-cell', 'A/Ohio/1/2013-egg', 'A/Romania/208058/2017-egg', 'A/SaoPaulo/76607/2015-egg', 'A/Kalamata/540/2017-cell', 'A/Kalamata/540/2017-egg', 'A/Mexico/InDRE2118/2005-egg', 'A/Ireland/10586/1999-egg', 'A/Shandong-Laicheng/1760/2016-egg', 'A/Qingdao/505/2009-egg', 'A/Corsica/3618-07/2015-egg', 'A/Missouri/1e/2005-egg', 'A/Parana/492/2017-egg', 'A/Bilthoven/908/1969', 'A/Bilthoven/908/1969-egg', 'A/CzechRepublic/46/2015-egg', 'A/Victoria/503/2015-egg', 'A/Victoria/503/2014-cell', 'A/Victoria/503/2014-egg', 'A/HongKong/1-6-MA21-1/1968', 'A/HongKong/1-6-MA21-1/1968-egg', 'A/Mexico/InDRE2160/2005-egg', 'A/Saitama/103/2014', 'A/Saitama/103/2014-cell', 'A/Saitama/103/2014-egg', 'A/SaoPaulo/13-13137/2013-egg', 'A/SaoPaulo/537178/2017-egg', 'A/SaoPaulo/15169/2014-egg', 'A/HoChiMinh/7218/2010-egg', 'A/Brisbane/190/2017-cell', 'A/Brisbane/190/2017-egg', 'A/HongKong/4801/2014-cell', 'A/HongKong/4801/2014-egg', 'A/Colombo/166/2003-egg', 'A/SaoPaulo/11-13826/2014-egg', 'A/Singapore/H2010370C/2010-egg', 'A/Kentucky/13e/2005-egg', 'A/Wyoming/9e/2005-egg', 'A/SaoPaulo/11194/2014-egg', 'A/Singapore/H2011140C/2011-egg', 'A/HoChiMinh/6318/2010-egg', 'A/SaoPaulo/55031/2014-egg', 'A/Tocantins/4622/2015-egg', 'A/CzechRepublic/38/2015-egg', 'A/Townsville/87/2010', 'A/Townsville/87/2010-cell', 'A/Townsville/87/2010-egg', 'A/HoChiMinh/5527/2010-egg', 'A/SaoPaulo/591445/2017-egg', 'A/Singapore/C2011573/2011-cell', 'A/Singapore/C2011573/2011-egg', 'A/Romania/802/2003-egg', 'A/SaoPaulo/9879/2014-egg', 'A/Singapore/H2011808bC/2011-egg', 'A/Corsica/3490-05/2015-egg', 'A/Jiangsu-Xinpu/1385/2012-egg', 'A/CzechRepublic/43/2015-egg', 'A/Romania/972/2003-egg', 'A/Colorado/6/2009-egg', 'A/SantaCatarina/451/2017-egg', 'A/MatoGrossodoSul/85939/2015-egg', 'A/Victoria/3/1999-cell', 'A/Victoria/3/1999-egg', 'A/Sydney/142/2016-cell', 'A/Sydney/142/2016-egg', 'A/Romania/3762-3764/2016-egg', 'A/CzechRepublic/41/2015-egg', 'A/Romania/206973/2017-egg', 'A/Hawaii/25/2012-cell', 'A/Hawaii/25/2012-egg', 'A/Goias/16075/2014-egg', 'A/Catalonia/NSVH110379682/2016-egg', 'A/Utah/7/2013-cell', 'A/Utah/7/2013-egg', 'A/RioGrandedoSul/671/2017-egg', 'A/Singapore/C2011803/2011-cell', 'A/Singapore/C2011803/2011-egg', 'A/Jiangsu-Sucheng/1309/2016-egg', 'A/Perth/98/2013-cell', 'A/Perth/98/2013-egg', 'A/Catalonia/NSVH100483235/2016-egg', 'A/MatoGrossodoSul/69756/2014-egg', 'A/Singapore/H2010321C/2010-egg', 'A/Jiangsu-Chongchuan/12179/2014-egg', 'A/Goias/852896/2018-egg', 'A/HongKong/1-12/1968-egg', 'A/SaoPaulo/52831/2014-egg', 'A/PortChalmers/JY2/1973-egg', 'A/Guangdong-Luohu/1257/2009', 'A/Guangdong-Luohu/1257/2009-egg', 'A/Ontario/C653602/2011-egg', 'A/SaoPaulo/419512/2017-egg', 'A/Wellington/2/2008', 'A/Wellington/2/2008-cell', 'A/Wellington/2/2008-egg', 'A/SaoPaulo/35339/2015-egg', 'A/Victoria/653/2017', 'A/Victoria/653/2017-cell', 'A/Victoria/653/2017-egg', 'A/Singapore/H2013422d/2013-egg', 'A/SaoPaulo/44541/2014-egg', 'A/Brisbane/32/2017-cell', 'A/Brisbane/32/2017-egg', 'A/Taiwan/1064/1985-egg', 'A/Jiangsu-Xinpu/1963/2014-egg', 'A/Singapore/H2013238/2013-egg', 'A/Georgia/4e/2006-egg', 'A/Taiwan/426/2012-egg', 'A/Nevada/20/2012-cell', 'A/Nevada/20/2012-egg', 'A/SouthAfrica/21/1998-egg', 'A/Newcastle/1/2013-cell', 'A/Newcastle/1/2013-egg', 'A/Ontario/C575475/2010-egg', 'A/SaoPaulo/58783/2014-egg', 'A/MatoGrosso/65925/2015-egg', 'A/HongKong/1-1-MA-20A/1968-egg', 'A/EspiritoSanto/58/2017-egg', 'A/Singapore/C2009485b/2009-cell', 'A/Singapore/C2009485b/2009-egg', 'A/Thessaloniki/62/2018-egg', 'A/Memphis/1/1968-egg', 'A/Indiana/52/2012-cell', 'A/Indiana/52/2012-egg', 'A/HoChiMinh/5677/2010-egg', 'A/RioGrandedoSul/663/2017-egg', 'A/Pennsylvania/CHOP4/2013-cell', 'A/Pennsylvania/CHOP4/2013-egg', 'A/SaoPaulo/496292/2017-egg', 'A/SouthAustralia/209/2017-cell', 'A/SouthAustralia/209/2017-egg', 'A/SaoPaulo/3-26596/2014-egg', 'A/Victoria/577/1999-egg', 'A/NewYork/39e/2004-egg', 'A/SaoPaulo/2-25786/2014-egg', 'A/Wyoming/3/2003', 'A/Wyoming/3/2003-egg', 'A/Victoria/504/2013-cell', 'A/Victoria/504/2013-egg', 'A/Singapore/C2011477/2011-cell', 'A/Singapore/C2011477/2011-egg', 'A/Victoria/507/2004', 'A/Victoria/507/2004-egg', 'A/Ohio/18/2012-cell', 'A/Ohio/18/2012-egg', 'A/SantaCatarina/455/2017-cell', 'A/SantaCatarina/455/2017-egg', 'A/Ontario/C537546/2010-egg', 'A/Kanagawa/IC1710/2017-cell', 'A/Kanagawa/IC1710/2017-egg', 'A/SaoPaulo/12-11047/2014-egg', 'A/Romania/206492/2017-egg', 'A/SouthAfrica/6388/2013-egg', 'A/Singapore/215/2008', 'A/Singapore/215/2008-egg', 'A/Taiwan/70307/2007-egg', 'A/SaoPaulo/11-28129/2014-egg', 'A/SaoPaulo/54152/2014-egg', 'A/Brisbane/321/2016-cell', 'A/Brisbane/321/2016-egg', 'A/SaoPaulo/428610/2017-egg', 'A/Anhui/1247/2005', 'A/Anhui/1247/2005-egg', 'A/Corsica/3500-10/2015-egg', 'A/HongKong/1-6/1968-egg', 'A/Singapore/H2011751/2011-cell', 'A/Singapore/H2011751/2011-egg', 'A/Brisbane/1/2013-cell', 'A/Brisbane/1/2013-egg', 'A/Brisbane/1/2012-cell', 'A/Brisbane/1/2012-egg', 'A/Brisbane/1/2017-cell', 'A/Brisbane/1/2017-egg', 'A/Brisbane/1/2018-cell', 'A/Brisbane/1/2018-egg', 'A/Virginia/5/2012-cell', 'A/Virginia/5/2012-egg', 'A/SaoPaulo/23630/2014-egg', 'A/Mexico/InDRE940/2003-egg', 'A/Brisbane/4/2007', 'A/Brisbane/4/2007-egg', 'A/Taiwan/240/2009-egg', 'A/Ontario/C582096/2010-egg', 'A/Perth/2/2014-egg', 'A/Virginia/16/2012-cell', 'A/Virginia/16/2012-egg', 'A/AthensGR/59/2012-cell', 'A/AthensGR/59/2012-egg', 'A/Taiwan/502/2013-egg', 'A/SaoPaulo/12-8700/2014-egg', 'A/CzechRepublic/64/2015-egg', 'A/SaoPaulo/53998/2014-egg', 'A/Corsica/3501-08/2015-egg', 'A/Montana/5/2011-cell', 'A/Montana/5/2011-egg', 'A/Ontario/R3296/2011-egg', 'A/HongKong/1-10-MA21-1/1968', 'A/HongKong/1-10-MA21-1/1968-egg', 'A/Georgia/168/2017-egg', 'A/Pennsylvania/CHOP3/2013-cell', 'A/Pennsylvania/CHOP3/2013-egg', 'A/SaoPaulo/5-1113/2014-egg', 'A/Singapore/C2011647/2011-cell', 'A/Singapore/C2011647/2011-egg', 'A/SaoPaulo/2-17089/2015-egg', 'A/SaoPaulo/6-15031/2015-egg', 'A/SaoPaulo/46667/2015-egg', 'A/SaoPaulo/5129/2014-egg', 'A/HoChiMinh/5186/2010-egg', 'A/HoChiMinh/7188/2010-egg', 'A/SaoPaulo/22896/2014-egg', 'A/England/1972-egg', 'A/Arizona/9/2012-cell', 'A/Arizona/9/2012-egg', 'A/HongKong/1-5-MA21-3/1968', 'A/HongKong/1-5-MA21-3/1968-egg', 'A/HongKong/1-2-MA21-3/1968', 'A/HongKong/1-2-MA21-3/1968-egg', 'A/Newcastle/22/2014-cell', 'A/Newcastle/22/2014-egg', 'A/SaoPaulo/72316/2015-egg', 'A/MatogrossodoSul/29492/2015-egg', 'A/SaoPaulo/23632/2014-egg', 'A/SaoPaulo/6-14445/2015-egg', 'A/Arkansas/18/2017', 'A/Arkansas/18/2017-egg', 'A/Taiwan/1091/1987-egg', 'A/HongKong/1-7/1968-egg', 'A/SaoPaulo/18819/2014-egg', 'A/Taiwan/70002/2007-egg', 'A/Shaanxi-Hanbin/11305/2016-cell', 'A/Shaanxi-Hanbin/11305/2016-egg', 'A/Taiwan/70120/2008-egg', 'A/Sofia/319/2007', 'A/Sofia/319/2007-egg', 'A/Sofia/682/2004', 'A/Sofia/682/2004-egg', 'A/Brazil/1898/2006-egg', 'A/Iowa/8/2011', 'A/Iowa/8/2011-cell', 'A/Iowa/8/2011-egg', 'A/NewHampshire/16/2012-egg', 'A/SaoPaulo/3-36687/2015-egg', 'A/Anhui/1238/2005', 'A/Anhui/1238/2005-egg', 'A/Mexico/InDRE835/2003-egg', 'A/Goias/57975/2014-egg', 'A/Brisbane/54/2009', 'A/Brisbane/54/2009-egg', 'A/SaoPaulo/51530/2014-egg', 'A/SaoPaulo/10-23262/2015-egg', 'A/Norway/1292/2012', 'A/Norway/1292/2012-egg', 'A/HoChiMinh/5947/2010-egg', 'A/DistritoFederal/756373/2018-egg', 'A/Geneva/4149722/2010', 'A/Geneva/4149722/2010-egg', 'A/Victoria/209/2009', 'A/Victoria/209/2009-egg', 'A/Fujian-Tongan/196/2009', 'A/Fujian-Tongan/196/2009-cell', 'A/Fujian-Tongan/196/2009-egg', 'A/Georgia/714/2017-egg', 'A/HongKong/1-11-MA21-1/1968-egg', 'A/Brisbane/220/2010', 'A/Brisbane/220/2010-cell', 'A/Brisbane/220/2010-egg', 'A/Perth/15/2009', 'A/Perth/15/2009-egg', 'A/Ontario/C5828/2011-egg', 'A/SaoPaulo/618711/2017-egg', 'A/MatoGrossodoSul/69755/2014-egg', 'A/Nepal/921/2006', 'A/Nepal/921/2006-egg', 'A/Tocantins/448840/2017-egg', 'A/SaoPaulo/9-130/2014-egg', 'A/SaoPaulo/17028/2014-egg', 'A/Taiwan/72106/2007-egg', 'A/Brisbane/6/2012-cell', 'A/Brisbane/6/2012-egg', 'A/Singapore/H2010619/2010-cell', 'A/Singapore/H2010619/2010-egg', 'A/DistrictofColumbia/1/2004-egg', 'A/Tocantins/18288/2015-egg', 'A/MatoGrossodoSul/45230/2014-egg', 'A/Corsica/15-07/2015-egg', 'A/SouthAfrica/8324/2013-egg', 'A/Taiwan/521/1980-egg', 'A/Ontario/P21546/2010-egg', 'A/SaoPaulo/604729/2017-egg', 'A/Fiji/4/2015-cell', 'A/Fiji/4/2015-egg', 'A/Pennsylvania/8/2013-cell', 'A/Pennsylvania/8/2013-egg', 'A/Anhui/1239/2005', 'A/Anhui/1239/2005-egg', 'A/Washington/16/2017', 'A/Washington/16/2017-cell', 'A/Washington/16/2017-egg', 'A/SaoPaulo/3-32837/2015-egg', 'A/SaoPaulo/461770/2017-egg', 'A/Victoria/8/2010', 'A/Victoria/8/2010-cell', 'A/Victoria/8/2010-egg', 'A/SaoPaulo/94077/2013-egg', 'A/Malaysia/1/2004', 'A/Malaysia/1/2004-cell', 'A/Malaysia/1/2004-egg', 'A/Mexico/InDRE2227/2005-egg', 'A/Ontario/C575590/2010-egg', 'A/HongKong/1-1-MA-20E/1968-egg', 'A/HongKong/1-1-MA-20B/1968-egg', 'A/HongKong/1-2/1968-egg', 'A/Taiwan/1110/1989-egg', 'A/CzechRepublic/120/2012-egg', 'A/Singapore/C2011452/2011-cell', 'A/Singapore/C2011452/2011-egg', 'A/Alaska/22/2012-cell', 'A/Alaska/22/2012-egg', 'A/Corsica/3490-2/2015-egg', 'A/Tocantins/448825/2017-egg', 'A/SouthAustralia/90/2017-cell', 'A/SouthAustralia/90/2017-egg', 'A/NewCaledonia/71/2014-egg', 'A/Santiago/7981/2006-egg', 'A/Kanagawa/IC1608/2017', 'A/Kanagawa/IC1608/2017-cell', 'A/Kanagawa/IC1608/2017-egg', 'A/Taiwan/9/2004-egg', 'A/Brisbane/29/2017-cell', 'A/Brisbane/29/2017-egg', 'A/Brisbane/185/2017-cell', 'A/Brisbane/185/2017-egg', 'A/Wyoming/5e/2005-egg', 'A/Hanoi/N015/2007-egg', 'A/CzechRepublic/26/2015-egg', 'A/Arizona/17/2013-cell', 'A/Arizona/17/2013-egg', 'A/SaoPaulo/13-8445/2014-egg', 'A/Memphis/3/1986-egg', 'A/Jiangsu-Sucheng/1489/2017-egg', 'A/Victoria/362/2011-cell', 'A/Victoria/362/2011-egg', 'A/England/42/1972', 'A/England/42/1972-cell', 'A/England/42/1972-egg', 'A/NewMexico/9/2014-cell', 'A/NewMexico/9/2014-egg', 'A/Sydney/30/2014-egg', 'A/Singapore/33/2009', 'A/Singapore/33/2009-cell', 'A/Singapore/33/2009-egg', 'A/Hawaii/44/2017', 'A/Hawaii/44/2017-cell', 'A/Hawaii/44/2017-egg', 'A/Corsica/40-09/2015-egg', 'A/Kanagawa/AC1619/2017-cell', 'A/Kanagawa/AC1619/2017-egg', 'A/SaoPaulo/6-10266/2014-egg', 'A/Singapore/Infimh-16-0019/2016-cell', 'A/Singapore/Infimh-16-0019/2016-egg', 'A/Taiwan/10/1979-egg', 'A/Victoria/527/2012-egg', 'A/Hanoi/BM38/2007-egg', 'A/CzechRepublic/60/2015-egg', 'A/SaoPaulo/44383/2014-egg', 'A/Sydney/52/2017-cell', 'A/Sydney/52/2017-egg', 'A/Sofia/141/2003', 'A/Sofia/141/2003-egg', 'A/Corsica/4069-10/2015-egg', 'A/SaoPaulo/54146/2014-egg', 'A/Taiwan/546/1980-egg', 'A/Mexico/InDRE756/2003-egg', 'A/Florida/2e/2006-egg', 'A/SriLanka/93/2017-cell', 'A/SriLanka/93/2017-egg', 'A/Johannesburg/3/1998-egg', 'A/SaoPaulo/3-28772/2014-egg', 'A/SriLanka/56/2013-egg', 'A/SaoPaulo/50552/2014-egg', 'A/Romania/206974/2017-egg', 'A/Southaustralia/9/2015-egg', 'A/Ohio/20/2012-cell', 'A/Ohio/20/2012-egg', 'A/HongKong/1186/2003-egg', 'A/Taiwan/1024/1983-egg', 'A/Jilin/17/2003', 'A/Jilin/17/2003-egg', 'A/MatoGrosso/71252/2015-egg', 'A/SantaCatarina/453/2017-egg', 'A/HoChiMinh/6117/2010-egg', 'A/Rondonia/5266/2015-egg', 'A/Tucuman/I2377/2014-egg', 'A/Singapore/H2010384/2010-cell', 'A/Singapore/H2010384/2010-egg', 'A/Tasmania/219/2016', 'A/Tasmania/219/2016-cell', 'A/Tasmania/219/2016-egg', 'A/Parana/718/2017-egg', 'A/CzechRepublic/45/2015-egg', 'A/Pennsylvania/14/2010', 'A/Pennsylvania/14/2010-cell', 'A/Pennsylvania/14/2010-egg', 'A/Victoria/746/2017', 'A/Victoria/746/2017-egg', 'A/Taiwan/1054/1985-egg', 'A/Romania/206618/2017-egg', 'A/Fiji/2/2015-cell', 'A/Fiji/2/2015-egg', 'A/CzechRepublic/18/2014-egg', 'A/HongKong/1-1-MA21-3/1968-egg', 'A/HongKong/4547/2009', 'A/HongKong/4547/2009-cell', 'A/HongKong/4547/2009-egg', 'A/EspiritoSanto/503/2017-egg', 'A/Texas/50/2012-cell', 'A/Texas/50/2012-egg', 'A/MatoGrossodoSul/66636/2014-egg', 'A/Iowa/119/2010-egg', 'A/Zagreb/1216/2007', 'A/Zagreb/1216/2007-egg', 'A/Nanchang/933/1995', 'A/Nanchang/933/1995-cell', 'A/Nanchang/933/1995-egg', 'A/SaoPaulo/11562/2014-egg', 'A/Taiwan/936/1981-egg', 'A/Singapore/GP637/2013-cell', 'A/Singapore/GP637/2013-egg', 'A/HoChiMinh/6398/2010-egg', 'A/HongKong/1-10/1968-egg', 'A/Brazil/1742/2005', 'A/Brazil/1742/2005-egg', 'A/SaoPaulo/625769/2017-egg', 'A/SouthCarolina/24/2012-cell', 'A/SouthCarolina/24/2012-egg', 'A/Washington/18/2013-cell', 'A/Washington/18/2013-egg', 'A/SouthAustralia/2/2013-egg', 'A/SaoPaulo/12-11416/2015-egg', 'A/CzechRepublic/25/2015-cell', 'A/CzechRepublic/25/2015-egg', 'A/Hunan-Beihu/1143/2011-egg', 'A/SaoPaulo/3-35911/2015-egg', 'A/Piaui/494713/2017-egg', 'A/MatoGrossodoSul/41977/2014-egg', 'A/Singapore/H2010301/2010-cell', 'A/Singapore/H2010301/2010-egg', 'A/SaoPaulo/53873/2013-egg', 'A/SaoPaulo/61282/2015-egg', 'A/Perth/201/2008', 'A/Perth/201/2008-egg', 'A/SaoPaulo/669992/2017-egg', 'A/Texas/141/2003-egg', 'A/Corsica/19-01/2015-egg', 'A/Singapore/21/2008', 'A/Singapore/21/2008-egg', 'A/SaoPaulo/12-11656/2015-egg', 'A/Fujian/445/2003', 'A/Fujian/445/2003-egg', 'A/Leningrad/360/1986', 'A/Leningrad/360/1986-egg', 'A/HoChiMinh/6388/2010-egg', 'A/SouthAfrica/8314/2013-egg', 'A/Romania/206396/2017-egg', 'A/Singapore/H2009389b/2009-cell', 'A/Singapore/H2009389b/2009-egg', 'A/Singapore/H2009332C/2009-egg', 'A/SouthAfrica/8766/2013-egg', 'A/Brisbane/341/2014-cell', 'A/Brisbane/341/2014-egg', 'A/Chile/8299/2003', 'A/Chile/8299/2003-egg', 'A/SouthAfrica/57/2000-egg', 'A/SaoPaulo/22768/2014-egg', 'A/Wisconsin/68e/2005-egg', 'A/Corsica/4069-09/2015-egg', 'A/Ontario/C835942/2010-egg', 'A/CzechRepublic/27/2015-egg', 'A/Victoria/318/1998-egg', 'A/CzechRepublic/20/2014-cell', 'A/CzechRepublic/20/2014-egg', 'A/CzechRepublic/20/2012-egg', 'A/Taiwan/70113/2007-egg', 'A/Pennsylvania/9/2012-cell', 'A/Pennsylvania/9/2012-egg', 'A/Guangdong-Luohu/1256/2009', 'A/Guangdong-Luohu/1256/2009-egg', 'A/Victoria/500/2004', 'A/Victoria/500/2004-egg', 'A/SaoPaulo/551/2016-egg', 'A/Fujian/140/2000', 'A/Fujian/140/2000-egg', 'A/MatoGrossodoSul/523805/2017-egg', 'A/NewHampshire/10/2012-cell', 'A/NewHampshire/10/2012-egg', 'A/Taiwan/601/1980-egg', 'A/Victoria/1968-egg', 'A/HoChiMinh/4596/2010-egg', 'A/Denmark/47/2014-cell', 'A/Denmark/47/2014-egg', 'A/Singapore/C2011411/2011-cell', 'A/Singapore/C2011411/2011-egg', 'A/Singapore/C2009863/2009-cell', 'A/Singapore/C2009863/2009-egg', 'A/Corsica/3490-03/2015-egg', 'A/Romania/207393/2017-egg', 'A/Switzerland/7587197/2012-cell', 'A/Switzerland/7587197/2012-egg', 'A/SaoPaulo/81060/2014-egg', 'A/Victoria/505/2004', 'A/Victoria/505/2004-egg', 'A/SaoPaulo/10767/2014-egg', 'A/SaoPaulo/3-27018/2014-egg', 'A/Victoria/512/2005', 'A/Victoria/512/2005-egg', 'A/SouthAfrica/6627/2013-egg', 'A/Geneva/1514818/2009', 'A/Geneva/1514818/2009-egg', 'A/MinasGerais/699/2017-egg', 'A/Almaty/3277/2012-egg', 'A/Ontario/C575589/2010-egg', 'A/Tocantins/49746/2015-egg', 'A/Kanagawa/ZC1613/2017-cell', 'A/Kanagawa/ZC1613/2017-egg', 'A/Wyoming/7e/2005-egg', 'A/CzechRepublic/132/2012-egg', 'A/Romania/206216/2016-egg', 'A/Perth/5/2008', 'A/Perth/5/2008-egg', 'A/Kenya/1544/2008-cell', 'A/Kenya/1544/2008-egg', 'A/Singapore/C2011545C/2011-egg', 'A/HongKong/1-12-MA21-1/1968-egg', 'A/SaoPaulo/12-15454/2014-egg', 'A/Catalonia/576683S/2016-egg', 'A/Romania/208461/2017-egg', 'A/Chongqing/Yuzhong159/2007-egg', 'A/Brisbane/198/2013-egg', 'A/Taiwan/851/2010', 'A/Taiwan/851/2010-egg', 'A/Ontario/N91737/2010-egg', 'A/Thessaloniki/63/2013-egg', 'A/Sydney/602/2009-egg', 'A/HongKong/32/2013-cell', 'A/HongKong/32/2013-egg', 'A/Tasmania/1012/2015-cell', 'A/Tasmania/1012/2015-egg', 'A/Maryland/2/2012-cell', 'A/Maryland/2/2012-egg', 'A/HubeiWujiagang/1121/2006', 'A/HubeiWujiagang/1121/2006-egg', 'A/UnitedKingdom/1861/2003-egg', 'A/Sydney/5/1997', 'A/Sydney/5/1997-egg', 'A/HongKong/2154/2009-egg', 'A/Brisbane/53/2009', 'A/Brisbane/53/2009-egg', 'A/Brisbane/10/2008', 'A/Brisbane/10/2008-egg', 'A/SaoPaulo/3-18529/2015-egg', 'A/HongKong/1-9/1968-egg', 'A/Ontario/C687408/2010-egg', 'A/Bangkok/1/1979', 'A/Bangkok/1/1979-cell', 'A/Bangkok/1/1979-egg', 'A/Goias/383452/2017-egg', 'A/Wyoming/5/2013-cell', 'A/Wyoming/5/2013-egg', 'A/Victoria/511/2004', 'A/Victoria/511/2004-cell', 'A/Victoria/511/2004-egg', 'A/SaoPaulo/734612/2018-egg', 'A/Rondonia/9097/2015-egg', 'A/Christchurch/516/2017-cell', 'A/Christchurch/516/2017-egg', 'A/SaoPaulo/61715/2013-egg', 'A/Jiangsu-Sucheng/1834/2017-egg', 'A/SaoPaulo/10774/2014-egg', 'A/Minnesota/11/2010', 'A/Minnesota/11/2010-egg', 'A/NewHampshire/3/2006-cell', 'A/NewHampshire/3/2006-egg', 'A/Tocantins/16261/2015-egg', 'A/SaoPaulo/4343/2014-egg', 'A/Taiwan/1097/1987-egg', 'A/HongKong/1-10-MA21-2/1968', 'A/HongKong/1-10-MA21-2/1968-egg', 'A/Chile/8266/2003', 'A/Chile/8266/2003-egg', 'A/Corsica/3500-18/2015-egg', 'A/SouthAustralia/21/2015-cell', 'A/SouthAustralia/21/2015-egg', 'A/Alaska/11e/2004-egg', 'A/SouthAfrica/9061/2013-egg', 'A/Romania/3784-3835/2016-egg', 'A/Singapore/H2011460a/2011-cell', 'A/Singapore/H2011460a/2011-egg', 'A/Texas/140/2003', 'A/Texas/140/2003-egg', 'A/SaoPaulo/49056/2015-egg', 'A/RhodeIsland/1/2010', 'A/RhodeIsland/1/2010-cell', 'A/RhodeIsland/1/2010-egg', 'A/Victoria/3/1975', 'A/Victoria/3/1975-egg', 'A/Uruguay/34/2005-egg', 'A/NewHampshire/9/2004-egg', 'A/NewJersey/3e/2006-egg', 'A/MatoGrossodoSul/510074/2018-egg', 'A/Taiwan/930/1981-egg', 'A/SaoPaulo/23803/2014-egg', 'A/SaoPaulo/60623/2014-egg', 'A/Mexico/InDRE2554/2011-egg', 'A/MatoGrosso/52531/2015-egg', 'A/Romania/207465/2017-egg', 'A/Cairo/51/2012-egg', 'A/PuertoRico/38/2011-egg', 'A/Auckland/19/1996-cell', 'A/Auckland/19/1996-egg', 'A/SaoPaulo/53699/2013-egg', 'A/Victoria/5060/2014-cell', 'A/Victoria/5060/2014-egg', 'A/Romania/3762-3765/2016-egg', 'A/Singapore/H2013729/2013-egg', 'A/SouthAfrica/8/1998-egg', 'A/Taiwan/1151/1992-egg', 'A/CzechRepublic/32/2015-egg', 'A/SouthAfrica/1147/1996-egg', 'A/SaoPaulo/3-26852/2014-egg', 'A/Canberra/1/1996-egg', 'A/Singapore/C2011301/2011-cell', 'A/Singapore/C2011301/2011-egg', 'A/Goias/42035/2014-egg', 'A/Wyoming/3e5/2003-egg', 'A/Brisbane/192/2017-cell', 'A/Brisbane/192/2017-egg', 'A/CzechRepublic/130/2012-egg', 'A/Romania/205950/2016-egg', 'A/MatoGrossodoSul/51672/2015-egg', 'A/Brisbane/174/2015-egg', 'A/NewYork/55e/2004-egg', 'A/MatoGrossodoSul/51673/2015-egg', 'A/HongKong/1-1-MA21-1/1968-egg', 'A/DominicanRepublic/3668/2009', 'A/DominicanRepublic/3668/2009-cell', 'A/DominicanRepublic/3668/2009-egg', 'A/Fujian/411/2002', 'A/Fujian/411/2002-cell', 'A/Fujian/411/2002-egg', 'A/Memphis/1/1986-egg', 'A/Johannesburg/33/1994', 'A/Johannesburg/33/1994-egg', 'A/Austria/64/2011-egg', 'A/Ontario/C432471/2010-egg', 'A/SaoPaulo/24361/2014-egg', 'A/Goias/886029/2018-egg', 'A/Victoria/265/2014-cell', 'A/Victoria/265/2014-egg', 'A/MatoGrossodoSul/554594/2017-egg', 'A/CzechRepublic/37/2015-egg', 'A/SaoPaulo/20448/2015-egg', 'A/HongKong/1-9-MA21-2/1968-egg', 'A/Goias/58899/2015-egg', 'A/Thessaloniki/19/2009', 'A/Thessaloniki/19/2009-egg', 'A/Tocantins/29517/2015-egg', 'A/HongKong/1-1/1968-egg', 'A/Brisbane/27/2017-cell', 'A/Brisbane/27/2017-egg', 'A/SaoPaulo/3-35046/2015-egg', 'A/Victoria/624/2017', 'A/Victoria/624/2017-cell', 'A/Victoria/624/2017-egg', 'A/Christchurch/512/2012-egg', 'A/MatoGrossodoSul/53099/2013-egg', 'A/SouthAustralia/9/2015-cell', 'A/SouthAustralia/9/2015-egg', 'A/Pennsylvania/4/2007', 'A/Pennsylvania/4/2007-cell', 'A/Pennsylvania/4/2007-egg', 'A/SaoPaulo/3-26965/2014-egg', 'A/SaoPaulo/7-15273/2013-egg', 'A/Auckland/21/2013-cell', 'A/Auckland/21/2013-egg', 'A/HongKong/1-1-MA21-2/1968-egg', 'A/CzechRepublic/121/2012-egg', 'A/Singapore/H2011507/2011-cell', 'A/Singapore/H2011507/2011-egg', 'A/Victoria/186/1982-egg', 'A/SaoPaulo/22930/2014-egg', 'A/SaoPaulo/3-29582/2015-egg', 'A/Singapore/H2013721a/2013-egg', 'A/Corsica/23-01/2015-egg', 'A/Corsica/4076-15/2015-egg', 'A/SaoPaulo/554097/2017-egg', 'A/SaoPaulo/21685/2015-egg', 'A/Victoria/700/2013', 'A/Victoria/700/2013-egg', 'A/Ontario/C728447/2010-egg', 'A/CzechRepublic/114/2012-egg', 'A/Hawaii/26/2012-cell', 'A/Hawaii/26/2012-egg', 'A/Nebraska/1/2018', 'A/Nebraska/1/2018-egg', 'A/SaoPaulo/35816/2014-egg', 'A/SaoPaulo/11199/2014-egg', 'A/Jiangsu-Qingpu/1660/2015-egg', 'A/SaoPaulo/10775/2014-egg', 'A/Corsica/N3849-04/2015-egg', 'A/Tasmania/37/2017-cell', 'A/Tasmania/37/2017-egg', 'A/Mexico/InDRE2246/2005-egg', 'A/Goias/85156/2014-egg', 'A/Kanagawa/IC1734/2018-cell', 'A/Kanagawa/IC1734/2018-egg', 'A/Taiwan/1127/1990-egg', 'A/Corsica/4076-12/2015-egg', 'A/HongKong/1-4-MA21-2/1968-egg', 'A/Taiwan/1452/2007-egg', 'A/DistritoFederal/21817/2014-egg', 'A/Jiangsu-Tianning/1707/2013-egg', 'A/Delaware/32/2016', 'A/Delaware/32/2016-cell', 'A/Delaware/32/2016-egg', 'A/SaoPaulo/446890/2017-egg', 'A/Hanoi/TX023/2007-egg', 'A/HongKong/2210/2009-cell', 'A/HongKong/2210/2009-egg', 'A/VIctoria/673/2014-egg', 'A/HongKong/1-9-MA21-3/1968-egg', 'A/Wisconsin/4/2018', 'A/Wisconsin/4/2018-cell', 'A/Wisconsin/4/2018-egg', 'A/Arkansas/13/2013-cell', 'A/Arkansas/13/2013-egg', 'A/SouthAustralia/135/2016-cell', 'A/SouthAustralia/135/2016-egg', 'A/Perth/16/2009', 'A/Perth/16/2009-cell', 'A/Perth/16/2009-egg', 'A/Brisbane/21/2013-cell', 'A/Brisbane/21/2013-egg', 'A/SaoPaulo/95739/2014-egg', 'A/Christchurch/28/2011-cell', 'A/Christchurch/28/2011-egg', 'A/Victoria/1000/2015-cell', 'A/Victoria/1000/2015-egg', 'A/SouthAfrica/7279/2013-egg', 'A/Kanagawa/AC1627/2017-cell', 'A/Kanagawa/AC1627/2017-egg', 'A/Madagascar/153/2010-egg', 'A/Florida/68/2012-egg', 'A/Alaska/24/2012-cell', 'A/Alaska/24/2012-egg', 'A/Romania/206437/2017-egg', 'A/Singapore/H2013721f/2013-egg', 'A/HongKong/JY2/1968-egg', 'A/Parana/487/2017-egg', 'A/Sichuan/346/1998-egg', 'A/Corsica/4041-02/2015-egg', 'A/Ontario/K4022160/2010-egg', 'A/Oregon/1/2013-cell', 'A/Oregon/1/2013-egg', 'A/CzechRepublic/22/2014-egg', 'A/Pennsylvania/CHOP8/2013-cell', 'A/Pennsylvania/CHOP8/2013-egg', 'A/Ontario/C76206/2011-egg', 'A/Singapore/H2011518/2011-cell', 'A/Singapore/H2011518/2011-egg', 'A/HongKong/1-2-MA21-1/1968-egg', 'A/CzechRepublic/93/2015-egg', 'A/Qingdao/522/2009-egg', 'A/HongKong/1-11-MA21-2/1968-egg', 'A/Singapore/H2011808a/2011-cell', 'A/Singapore/H2011808a/2011-egg', 'A/Colombo/143/2003-egg', 'A/Corsica/4069-18/2015-egg', 'A/DistritoFederal/524550/2017-egg', 'A/Goias/82117/2014-egg', 'A/SaoPaulo/505950/2017-egg', 'A/Ohio/17/2012-cell', 'A/Ohio/17/2012-egg', 'A/Auckland/5/1996', 'A/Auckland/5/1996-cell', 'A/Auckland/5/1996-egg', 'A/SaoPaulo/28596/2015-egg', 'A/SaoPaulo/913733/2017-egg', 'A/Wyoming/28e/2005-egg', 'A/Ohio/23/2012-cell', 'A/Ohio/23/2012-egg', 'A/SaoPaulo/60985/2015-egg', 'A/SaoPaulo/6-10330/2014-egg', 'A/Quebec/26-191206/2006-egg', 'A/Kentucky/3e/2006-egg', 'A/Pennsylvania/CHOP2/2013-cell', 'A/Pennsylvania/CHOP2/2013-egg', 'A/Netherlands/938/1992', 'A/Netherlands/938/1992-egg', 'A/Hanoi/TX058/2007-egg', 'A/SaoPaulo/10-20620/2015-egg', 'A/Ontario/C543493/2010-egg', 'A/Ontario/C575591/2010-egg', 'A/Brisbane/285/2016-cell', 'A/Brisbane/285/2016-egg', 'A/Taiwan/1149/1992-egg', 'A/Pennsylvania/15/2010', 'A/Pennsylvania/15/2010-cell', 'A/Pennsylvania/15/2010-egg', 'A/Cardiff/7109/2018-egg', 'A/Norway/3806/2016', 'A/Norway/3806/2016-cell', 'A/Norway/3806/2016-egg', 'A/Singapore/GP1940/2013-cell', 'A/Singapore/GP1940/2013-egg', 'A/Sichuan/262/2003-egg', 'A/Washington/106/2016', 'A/Washington/106/2016-cell', 'A/Washington/106/2016-egg', 'A/Victoria/563/2010', 'A/Victoria/563/2010-cell', 'A/Victoria/563/2010-egg', 'A/Ontario/C526567/2010-egg', 'A/SaoPaulo/10-38324/2013-egg', 'A/Tehran/83/1999-egg', 'A/SouthAustralia/30/2012-egg', 'A/Washington/17/2013-cell', 'A/Washington/17/2013-egg', 'A/Goias/346863/2017-egg', 'A/SaoPaulo/3-34891/2014-egg', 'A/SaoPaulo/72317/2015-egg', 'A/HongKong/1-11-MA21-3/1968-egg', 'A/SaoPaulo/7-6864/2014-egg', 'A/Kanagawa/ZC1709/2017-cell', 'A/Kanagawa/ZC1709/2017-egg', 'A/Indiana/10/2011-cell', 'A/Indiana/10/2011-egg', 'A/Wyoming/1e/2006-egg', 'A/Iowa/19/2010-cell', 'A/Iowa/19/2010-egg', 'A/EspiritoSanto/506/2017-egg', 'A/HongKong/1-1-MA-12/1968-egg', 'A/SouthAustralia/273/2017-cell', 'A/SouthAustralia/273/2017-egg', 'A/Pennsylvania/CHOP5/2013-egg', 'A/Victoria/523/2004', 'A/Victoria/523/2004-cell', 'A/Victoria/523/2004-egg', 'A/SaoPaulo/55257/2014-egg', 'A/SaoPaulo/48734/2014-egg', 'A/Palau/6720/2014-cell', 'A/Palau/6720/2014-egg', 'A/CzechRepublic/98/2012-egg', 'A/Goias/57978/2014-egg', 'A/SaoPaulo/691063/2018-egg', 'A/Perth/501/2010', 'A/Perth/501/2010-cell', 'A/Perth/501/2010-egg', 'A/Pennsylvania/CHOP7/2013-cell', 'A/Pennsylvania/CHOP7/2013-egg', 'A/Canberra/82/2014-cell', 'A/Canberra/82/2014-egg', 'A/PortChalmers/1/1973-cell', 'A/PortChalmers/1/1973-egg', 'A/SaoPaulo/53849/2014-egg', 'A/HongKong/1-2-MA21-2/1968', 'A/HongKong/1-2-MA21-2/1968-egg', 'A/Jiangsu-Chongchuan/12714/2015-egg', 'A/SaoPaulo/2-15508/2015-egg', 'A/AmericanSamoa/4786/2013-cell', 'A/AmericanSamoa/4786/2013-egg', 'A/Singapore/H2009471C/2009-egg', 'A/SaoPaulo/690385/2018-egg', 'A/PuertoRico/36/2011-cell', 'A/PuertoRico/36/2011-egg', 'A/PuertoRico/36/2012-cell', 'A/PuertoRico/36/2012-egg', 'A/Kanagawa/AC1736/2018-egg', 'A/SaoPaulo/65194/2015-egg', 'A/SaoPaulo/596258/2017-egg', 'A/Fes/314/1998-egg', 'A/SaoPaulo/49665/2015-egg', 'A/SaoPaulo/72626/2014-egg', 'A/Delaware/23/2012-cell', 'A/Delaware/23/2012-egg', 'A/Mexico/InDRE2601/2003-egg', 'A/Chita/6/2003', 'A/Chita/6/2003-egg', 'A/Georgia/278/2017-egg', 'A/Serres/77/2007', 'A/Serres/77/2007-egg', 'A/Kenya/118/2013-cell', 'A/Kenya/118/2013-egg', 'A/SaoPaulo/50798/2015-egg', 'A/SaoPaulo/1044167/2017-egg', 'A/CzechRepublic/19/2014-cell', 'A/CzechRepublic/19/2014-egg', 'A/Newcastle/36/2017-cell', 'A/Newcastle/36/2017-egg', 'A/Norway/2178/2014', 'A/Norway/2178/2014-cell', 'A/Norway/2178/2014-egg', 'A/Georgia/3/2013-cell', 'A/Georgia/3/2013-egg', 'A/Tocantins/444385/2017-egg', 'A/Venezuela/6971/2005-egg', 'A/Shandong-Laicheng/1104/2016-egg', 'A/RioGrandedoSul/667/2017-egg', 'A/Singapore/H2009334C/2009-egg', 'A/MatoGrosso/44112/2015-egg', 'A/Hawaii/27/2012-cell', 'A/Hawaii/27/2012-egg', 'A/Victoria/210/2009', 'A/Victoria/210/2009-cell', 'A/Victoria/210/2009-egg', 'A/SaoPaulo/11-10263/2014-egg', 'A/Nebraska/19/2015', 'A/Nebraska/19/2015-cell', 'A/Nebraska/19/2015-egg', 'A/SantaCatarina/141/2017', 'A/SantaCatarina/141/2017-egg', 'A/Singapore/H2011797/2011-cell', 'A/Singapore/H2011797/2011-egg', 'A/Brisbane/11/2010', 'A/Brisbane/11/2010-cell', 'A/Brisbane/11/2010-egg', 'A/Brisbane/10/2007', 'A/Brisbane/10/2007-cell', 'A/Brisbane/10/2007-egg', 'A/SaoPaulo/3-29607/2015-egg', 'A/Taiwan/1069/1985-egg', 'A/Goias/44075/2015-egg', 'A/Singapore/C2011422/2011-cell', 'A/Singapore/C2011422/2011-egg', 'A/Goias/51166/2014-egg', 'A/CzechRepublic/35/2015-egg', 'A/MatoGrossodoSul/48293/2014-egg', 'A/Philippines/2-MA/1982-egg', 'A/SouthAustralia/54/2016-egg', 'A/Thailand/5107/2004-egg', 'A/Shanghai/369/2003', 'A/Shanghai/369/2003-egg', 'A/Vermont/28/2012-cell', 'A/Vermont/28/2012-egg', 'A/HongKong/1-9-MA21-1/1968-egg', 'A/SaoPaulo/40792/2015-egg', 'A/Brisbane/3e/2005-egg', 'A/Jilin/16/2003', 'A/Jilin/16/2003-egg', 'A/Singapore/C2011614/2011-cell', 'A/Singapore/C2011614/2011-egg', 'A/SaoPaulo/9-1757/2015-egg', 'A/Goias/67636/2014-egg', 'A/Romania/207297/2017-egg', 'A/Goias/49695/2015-egg', 'A/Guadeloupe/4127/2014-cell', 'A/Guadeloupe/4127/2014-egg', 'A/CzechRepublic/63/2015-egg', 'A/Mexico/InDRE2662/2003-egg', 'A/Wyoming/3e8/2003-egg', 'A/Wyoming/2/2006-egg', 'A/SaoPaulo/10-19619/2015-egg', 'A/HongKong/1-5-MA21-2/1968-egg', 'A/Taiwan/471/2008', 'A/Taiwan/471/2008-egg', 'A/Singapore/KK943/2014-egg', 'A/Goias/16069/2014-egg', 'A/Singapore/H2011471b/2011-cell', 'A/Singapore/H2011471b/2011-egg', 'A/Bangladesh/9334/2013-cell', 'A/Bangladesh/9334/2013-egg', 'A/Corsica/33-02/2015-egg', 'A/Kanagawa/ZC1617/2017-cell', 'A/Kanagawa/ZC1617/2017-egg', 'A/Parana/723/2017-egg', 'A/SaoPaulo/19850/2015-egg', 'A/Brazil/1499/2005', 'A/Brazil/1499/2005-egg', 'A/Texas/6/2014-cell', 'A/Texas/6/2014-egg', 'A/Corsica/4076-09/2015-egg', 'A/Bangladesh/5014/2009', 'A/Bangladesh/5014/2009-egg', 'A/Newcastle/35/2017-cell', 'A/Newcastle/35/2017-egg', 'A/HongKong/1-12-MA21-2/1968-egg', 'A/Taiwan/1143/1992-egg', 'A/Taiwan/929/1981-egg', 'A/Uruguay/716/2007', 'A/Uruguay/716/2007-cell', 'A/Uruguay/716/2007-egg', 'A/StPetersburg/5/2009-cell', 'A/StPetersburg/5/2009-egg', 'A/Bangladesh/9269/2013-cell', 'A/Bangladesh/9269/2013-egg', 'A/SouthAfrica/8806/2013-egg', 'A/Singapore/C2011244/2011-cell', 'A/Singapore/C2011244/2011-egg', 'A/SaoPaulo/21804/2015-egg', 'A/SaoPaulo/6-15209/2015-egg', 'A/SouthAustralia/55/2014-cell', 'A/SouthAustralia/55/2014-egg', 'A/Rondonia/492910/2017-egg', 'A/Wyoming/8e/2005-egg', 'A/Brisbane/299/2016-cell', 'A/Brisbane/299/2016-egg', 'A/Brisbane/299/2011-cell', 'A/Brisbane/299/2011-egg', 'A/Singapore/36/2004', 'A/Singapore/36/2004-egg']

# vaccine choices copied from https://github.com/blab/nextflu/tree/de1c8323f75b0fbac9cf26451380d2758288290e/auspice/_includes feb 1 2018
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
        'A/Singapore/Infimh-16-0019/2016': "2017-09-28"
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
    '12y': {"tau": 0.0005, "time_window": 0.5}
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

lineage_to_glyc_mask = {
    "h3n2": "ha1_h3n2",
    "h1n1pdm": "ha1_globular_head_h1n1pdm"
}

clade_designations = {
    "h3n2":{
        "3b":    [('HA2',158,'N'), ('HA1',198,'S'), ('HA1',312,'S'), ('HA1',223,'I'),
            ('HA1',145,'S')],
        "3c":    [('HA1',45,'N'), ('HA1',48,'I'), ('nuc',456,'T'), ('HA1',198,'S'),
            ('HA1',312,'S'), ('HA1',223,'I')],
        "3c2":   [('HA2',160,'N'), ('HA1',145,'S'), ('nuc',693,'A'), ('nuc',1518,'G')],
        "3c3":   [('nuc',285,'T'), ('nuc',430,'G'), ('nuc',472,'G'), ('nuc',1296,'A')],
        "3c3.A": [('HA1',128,'A'), ('HA1',142,'G'), ('HA1',159,'S')],
        "3c2.A": [('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'),
            ('HA1',311,'H'), ('nuc',1491,'A'), ('nuc', 234, 'A')],
        "A1":[('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'), # formerly clade 3c2.a1
            ('HA1',311,'H'), ('nuc',1491,'A'), ('HA1',171,'K'),
            ('HA2',77,'V'), ('HA2',155,'E'), ('HA2',160,'N')],
        "3c2":   [('HA1',144,'N'), ('HA1',159,'F'), ('HA1',225,'N'),
            ('HA1',160,'T'), ('HA1',142,'R')],
        "3c3.B": [('HA1',83,'R'), ('HA1',261,'Q'), ('HA1',62,'K'),
            ('HA1',122,'D')],
        "A2": [('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'), # formerly clade 3
            ('HA1',311,'H'), ('nuc',1491,'A'), ('nuc', 234, 'A'),
            ('HA1',131,'K'), ('HA1',142,'K'), ('HA1',261,'Q')],
        "A2/re": [('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'), # formerly clade 3
            ('HA1',311,'H'), ('nuc',1491,'A'), ('nuc', 234, 'A'),
            ('HA1',131,'K'), ('HA1',142,'K'), ('HA1',261,'Q'),
            ('nuc',1689,'T')],
        "A3": [('HA1',159,'Y'), ('HA1',225,'D'), ('HA1',311,'H'), # formerly clade 2
            ('nuc',1491,'A'), ('nuc',234,'A'), ('HA1',121,'K'),
            ('HA1',144,'K')],
        "A4": [('HA1',159,'Y'), ('HA1',225,'D'), ('HA1',311,'H'), # formerly clade 1
            ('nuc',1491,'A'), ('nuc',234,'A'), ('HA1',53,'N'),
            ('HA1',144,'R'), ('HA1',171,'K'), ('HA1',192,'T'),
            ('HA1',197,'H')],
        "A1a": [('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'), # formerly clade 4
            ('HA1',311,'H'), ('nuc',1491,'A'), ('nuc', 234, 'G'),
            ('HA2',150,'E'), ('nuc',114,'T')],
        "A1b": [('HA1',159,'Y'), ('HA1',225,'D'), ('nuc',1491,'A'), # formerly clade 5
            ('nuc', 234, 'G'), ('HA1',92,'R'), ('HA1',311,'Q'),
            ('nuc',538,'C')],
        "A1b/135N": [('HA1',159,'Y'), ('HA1',225,'D'), ('nuc',1491,'A'),
            ('nuc', 234, 'G'), ('HA1',92,'R'), ('HA1',311,'Q'),
            ('nuc',538,'C'), ('nuc',81,'G'), ('nuc',453,'T')],
        "A1b/135K": [('HA1',159,'Y'), ('HA1',225,'D'), ('nuc',1491,'A'),
            ('nuc', 234, 'G'), ('HA1',92,'R'), ('HA1',311,'Q'),
            ('nuc',538,'C'), ('nuc',233,'G'), ('nuc',472,'G'),
            ('nuc',452,'A')]
    },
    "h1n1pdm":{
        '1': [('HA1', 125, 'N'), ('HA1', 134 ,'A'), ('HA1', 183, 'S'), ('HA1', 31,'N'), ('HA1', 216, 'I')],
        '2': [('HA1', 125, 'N'), ('HA1', 134 ,'A'), ('HA1', 183, 'S'), ('HA1', 31,'D')],
        '3': [('HA1', 134 ,'T'), ('HA1', 183, 'P')],
        '4': [('HA1', 125, 'D'), ('HA1', 134 ,'A'), ('HA1', 183, 'S')],
        '5': [('HA1', 87, 'N'), ('HA1', 205, 'K'), ('HA1', 216, 'V'), ('HA1', 149, 'L')],
        '6': [('HA1', 185,'T'),  ('HA1', 97, 'N'), ('HA1', 197, 'A')],
        '6c':[('HA1', 234,'I'),  ('HA1', 97, 'N'), ('HA1', 197, 'A'), ('HA1', 283,'E')],
        '6b':[('HA1', 163,'Q'),  ('HA1', 256, 'T'), ('HA1', 197, 'A'), ('HA1', 283,'E')],
        '7': [('HA1', 143,'G'),  ('HA1', 97, 'D'), ('HA1', 197, 'T')],
        '8': [('HA1', 186,'T'),  ('HA1', 272,'A')],
        '6b.1':[('HA1', 163,'Q'),  ('HA1', 256, 'T'), ('HA1', 197, 'A'), ('HA1', 283, 'E'), ('SigPep', 13, 'T'), ('HA1', 84, 'N'), ('HA1', 162, 'N')],
        '6b.2':[('HA1', 163,'Q'),  ('HA1', 256, 'T'), ('HA1', 197, 'A'), ('HA1', 283, 'E'), ('HA2', 164, 'G'), ('HA1', 152, 'T'), ('HA2', 174, 'E')],
        '164T':[('HA1', 84, 'N'), ('HA1', 162, 'N'), ('HA1', 74, 'R'), ('HA1', 164, 'T'), ('HA1', 295, 'V')]
    },
    "vic":{
        '1A': [('nuc',206,'G'), ('nuc',644,'C'), ('nuc',1340,'T'), ('nuc',1821,'T'),
            ('HA1',165,'K'), ('HA1',172,'P')],
        '1B': [('nuc',1034,'G'), ('nuc',1172,'G'), ('HA1',165,'K'), ('HA1',172,'P')],
        '117V': [('HA1', 75,'K'), ('HA1', 58, 'L'), ('HA1', 165, 'K'), ('HA1', 129, 'D'), ('HA1', 117, 'V')],
        'DV': [('HA1', 75,'K'), ('HA1', 58, 'L'), ('HA1', 165, 'K'), ('HA1', 117, 'V'), ('HA1', 180, 'V'), ('HA2', 152, 'K'), ('HA1', 162, 'X')]
    },
    "yam":{
        '2':  [('HA1', 48,'K'), ('HA1', 108, 'A'), ('HA1', 150, 'S')],
        '3':  [('HA1', 48,'R'), ('HA1', 108, 'P'), ('HA1', 150, 'I')],
        '3a': [('HA1', 37,'A'), ('HA1', 298, 'E'), ('HA1', 48,'R'), ('HA1', 105, 'P'), ('HA1', 150, 'I')],
        '172Q': [('HA1', 48,'R'), ('HA1', 108, 'P'), ('HA1', 150, 'I'), ('HA1', 116, 'K'), ('HA1', 172, 'Q'), ('HA1', 298, 'E'), ('HA1', 312, 'K')]
    }
}
