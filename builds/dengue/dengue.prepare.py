from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import base.prepare
from base.prepare import prepare
from base.utils import fix_names
from base.titer_model import TiterCollection
from datetime import datetime, timedelta, date
import argparse
from dengue_subsampling import dengue_subsampling
import numpy as np

def collect_args():
    """Returns a dengue-specific argument parser.
    """
    parser = base.prepare.collect_args()

    parser.add_argument('-s', '--serotypes', '--lineage', choices=["all", "denv1", "denv2", "denv3", "denv4", "multiple"], default="multiple", type=str, nargs='+', help="Serotype(s) to prepare. \"multiple\" will run them all.")
    parser.add_argument('-y', '--years_back', type=int, default=100, help="Years back in time to sample from")
    parser.add_argument('--titers', default='../../../fauna/data/dengue_titers.tsv', help="tab-delimited file of titer strains and values from fauna (e.g., dengue_titers.tsv)")
    parser.add_argument('--strains', help="a text file containing a list of strains (one per line) to prepare without filtering or subsampling")
    parser.set_defaults(
        viruses_per_month=3,
        file_prefix=None
    )
    return parser

dropped_strains = [
    'DENV1/VIETNAM/BIDV992/2006', 'DENV1/FRANCE/00475/2008', 'DENV1/VIETNAM/BIDV3990/2008', 'DENV2/HAITI/DENGUEVIRUS2HOMOSAPIENS1/2016', # Probable recombinants
    'DENV2/AUSTRALIA/QML22/2015', # Suspiciously far diverged
    'DENV2/MALAYSIA/DKD811/2008', 'DENV2/MALAYSIA/P81407/1970', 'DENV2/SENEGAL/0674/1970', 'DENV2/SENEGAL/DAKAR0761/1974', # Sylvatic
    'DENV2/NIGERIA/IBH11234/1966', 'DENV2/NIGERIA/IBH11664/1966', 'DENV2/NIGERIA/IBH11208/1966', 'DENV2/SENEGAL/DAKARD75505/1999', # Sylvatic
    'DENV2/SENEGAL/DAKAR141069/1999', 'DENV2/SENEGAL/DAKAR141070/1999', 'DENV2/GUINEA/PM33974/1981', 'DENV2/BURKINA_FASO/DAKAR2039/1980', # Sylvatic
    'DENV2/COTE_D_IVOIRE/DAKAR578/1980', 'DENV2/COTE_D_IVOIRE/DAKAR510/1980', 'DENV2/MALAYSIA/SAB/2015', 'DENV2/TRINIDAD_AND_TOBAGO/NA/1953'# Sylvatic
    'DENV4/MALAYSIA/P731120/1973', 'DENV4/MALAYSIA/P215/1975'# Sylvatic
]

genotype_references = ['AB074760', 'EU848545', 'D10513', 'JF297570', 'JQ922547', 'JF297572', 'JQ922544', 'AF425629', 'AF180817', 'AF425625', 'JQ922546', 'AY713473', 'EF457905', 'U88535', 'AY722801', 'JN379478', 'AY732476', 'AY732474', 'AF425630', 'AY732483', 'AF425611', 'GQ868601', 'AF425620', 'AF425628', 'AB600923', 'AF226687', 'M87512', 'JN638342', 'AF226685', 'AY732477', 'AF425626', 'AY732413', 'FJ196845', 'DQ285562', 'AY153757', 'AY762084', 'AY732461', 'AB608789', 'AM746218', 'AY732480', 'FJ196846', 'JN638344', 'AF425632', 'JN379486', 'DQ855297', 'AY722802', 'JN638340', 'AF311956', 'AB189121', 'AF309641', 'AF298808', 'GQ868559', 'AY726555', 'AB189121', 'AF298807', 'AY422785', 'AF514883', 'AF514889', 'AY620951', 'JN415515', 'JN415499', 'DQ672564', 'AY713476', 'AY732479', 'AY732464', 'AY630407', 'DQ672564', 'KF559254', 'JN415503', 'EU863650', 'FJ196842', 'JN379472', 'JN415488', 'AB195673', 'EU596501', 'JN415524', 'FR666922', 'AM746216', 'AY835999', 'AY858983', 'AB204803', 'EU081258', 'EU081226', 'FN825674', 'EU179860', 'JF297581', 'JN415516', 'JQ922548', 'FJ196844', 'HQ332182', 'JN415533', 'JN415523', 'GU131962', 'KC172829', 'AB608787', 'JN415506', 'GU131863', 'JN903579', 'GQ357692', 'GU131792', 'JN415521', 'JN415527', 'JN415534', 'KR919821', 'KC692495', 'JN903581', 'JF960211', 'JN415519', 'GU131895', 'HQ891316', 'JQ915080', 'HG316481', 'HG316482', 'KC692512', 'KC182084', 'KR919811', 'JN415492', 'JN415494', 'JN415513', 'KJ189367', 'JQ675358', 'KJ649286', 'KR919805', 'KC848576', 'KR919815', 'JX298570', 'KF289072', 'KJ933413', 'KJ726662', 'KR919808', 'KR919819', 'KR919813', 'KF973455', 'KF864667', 'KR919806', 'KR919807', 'KR919809', 'KR919817', 'KF184975', 'KM458188', 'LC002828', 'KR919816', 'KR919812', 'KR919818', 'KR919814', 'KR919820', 'KR919810', 'AF038403', 'U87411', 'EF105387', 'EU003591', 'EF105388', 'EF105384', 'EF105379', 'DQ181806', 'EF105385', 'EU056812', 'EF105382', 'EF105380', 'EF105383', 'EF105386', 'EF105381', 'EF105378', 'M20558', 'AY037116', 'JX966379', 'DQ181799', 'AF208496', 'AF359579', 'EF105389', 'EF105390', 'EF457904', 'EU179858', 'FJ467493', 'JF260983', 'KC294207', 'M93130', 'AF317645', 'AY676353', 'JQ411814', 'AY744685', 'AY923865', 'DQ675519', 'JN406514', 'AY099337', 'AY099336', 'EF629367', 'EU081198', 'HQ332171', 'JN406515', 'HG316484', 'JF262783', 'JF262780', 'EF457906', 'JF262779', 'AF326573', 'JF262782', 'JF262781', 'AY618989', 'AY618993', 'AY618992', 'KF041260', 'JQ822247', 'JN983813', 'JX024758']

sanofi_vaccine_strains = {
    'denv1': 'DENV1/THAILAND/PUO359/1980',
    'denv2': 'DENV2/THAILAND/PUO218/1980',
    'denv3': 'DENV3/THAILAND/PAH88188/1988',
    'denv4': 'DENV4/INDONESIA/S1228/1978'}

references = {
    "denv1": {"metadata": {'strain': "DENV1/NAURUISLAND/REFERENCE/1997", "accession": "NC_001477", "date": "1997-XX-XX", 'host': "NA", 'country': "Nauru", 'region': "oceania"}},
    "denv2": {"metadata": {'strain': "DENV2/THAILAND/REFERENCE/1964", "accession": "NC_001474", "date": "1964-XX-XX", 'host': "NA", 'country': "Thailand", "region": "southeast_asia"}},
    "denv3": {"metadata": {'strain': "DENV3/SRI_LANKA/REFERENCE/2000", "accession": "NC_001475", "date": "2000-XX-XX", 'host': "NA", 'country': "Sri Lanka", "region": "south_asia"}},
    "denv4": {"metadata": {'strain': "DENV4/NA/REFERENCE/2003", "accession": "NC_002640", "date": "2003-XX-XX", 'host': "NA", 'country': "NA", "region": "NA"}},
}

for key in references:
    references[key]["genes"] = ['C', 'M', 'E', 'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', '2K', 'NS4B', 'NS5']
    references[key]["path"] = "metadata/dengue_{}_outgroup.gb".format(key)
    references[key]["include"] = 0
references["all"] = references["denv4"]


##### Utils #####
def select_serotype(infile, path, serotype):
    '''
    Takes all-serotype fasta file, path to save output, and desired serotype.
    Writes appropriate subset of sequences to dengue_serotype.fasta
    Returns path to output file as string.
    '''
    from Bio import SeqIO
    sequences = [ i for i in SeqIO.parse(infile, 'fasta') if i.description.split('/')[0] == serotype.upper() ]
    SeqIO.write(sequences, path+'dengue_%s.fasta'%serotype, 'fasta')
    return path+'dengue_%s.fasta'%serotype

##### Set up config #####
# For dengue, config is a function so it is applicable for multiple lineages
def make_config(serotype, params):
    if params.file_prefix is not None:
        file_prefix = params.file_prefix
    else:
        file_prefix = "dengue_%s" % serotype

    if params.sequences is not None:
        input_paths = [params.sequences]
    elif os.path.isfile("../../../fauna/data/dengue_%s.fasta"%serotype): #is file: # Look for a serotype-specific fasta
        input_paths = ["../../../fauna/data/dengue_%s.fasta"%serotype]
    else: # If it doesn't exist, try to pull serotype-specific sequences out of the all-serotype fasta (warn the user of this behavior)
        input_paths = [select_serotype('../../../fauna/data/dengue_all.fasta', '../../../fauna/data/', serotype)]
        print('WARNING: Did not find serotype-specific fasta file.\nPulled sequences with serotype %s from all-serotype fasta file %s\nWrote these to file %s'%(serotype, '../fauna/data/dengue.fasta', input_paths))

    years_back = params.years_back
    time_interval = [datetime.today().date(), (datetime.today()  - timedelta(days=365.25 * years_back)).date()]

    if params.titers is not None:
        if not os.path.isfile(params.titers):
            params.titers = '../../../fauna/data/%s'%params.titers
        titer_values, strains, sources = TiterCollection.load_from_file(params.titers)
    else:
        titer_values, strains, sources = None, None, None

    force_include = sanofi_vaccine_strains.values()

    config = {
        "dir": "dengue",
        "lineage": serotype,
        "title": "Real-time tracking of dengue evolution",
        "maintainer": ["Sidney Bell", "http://bedford.io/team/sidney-bell/"],
        "file_prefix": file_prefix,
        "input_paths": input_paths,
        "header_fields": {0:'strain', 1:'accession', 2:'date', 3:'region', 4:'country',
                        5:'division', 6: 'location', 7: 'authors', 8: 'url'},
        "filters": (("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
            ("Bad Region", lambda s: any([
                                        s.attributes['region'] not in ['', ' ', '?'],
                                        s.attributes['accession'] in force_include,
                                        s.attributes['strain'] in force_include
                                        ]))),

        "subsample": dengue_subsampling(params, years_back, titer_values,
        force_include),

        "add_urls": {
            "prefix": "https://www.ncbi.nlm.nih.gov/nuccore/%s",
            "attr": "accession" },

        "colors": ["authors", "region", "country"],
        "lat_longs": ["region", "country"],
        "auspice_filters": ["authors", "region", "country"],
        "reference": references[serotype],
        "time_interval": time_interval,
        "titers": titer_values,
        "strains": params.strains,
        "sources": sources
    }
    return config


if __name__=="__main__":
    parser = collect_args()
    params = parser.parse_args()

    if 'multiple' in params.serotypes:
        params.serotypes = ['denv1', 'denv2', 'denv3', 'denv4', 'all']

    # modify config file according to serotype
    for serotype in params.serotypes:
        print("Preparing serotype %s"%serotype)

        config = make_config(serotype, params)

        if params.viruses_per_month == 0:
            config["subsample"] = False
        else:
            config["subsample"]["threshold"] = params.viruses_per_month

        runner = prepare(config)
        runner.load_references()
        runner.applyFilters()
        runner.subsample()
        runner.colors()
        runner.latlongs()
        runner.write_to_json()
