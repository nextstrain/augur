from __future__ import print_function
import os, sys
sys.path.append('..')
import base.prepare
from base.prepare import prepare
from datetime import datetime
from base.utils import fix_names
from datetime import datetime, timedelta, date
import argparse

def collect_args():
    """Returns a dengue-specific argument parser.
    """
    parser = base.prepare.collect_args()

    parser.add_argument('-s', '--serotypes', '--lineage', choices=["all", "denv1", "denv2", "denv3", "denv4", "multiple"], default="multiple", type=str, nargs='+', help="Serotype(s) to prepare. \"multiple\" will run them all.")
    parser.add_argument('-y', '--years_back', type=int, default=100, help="Years back in time to sample from")
    parser.set_defaults(
        viruses_per_month=3,
        file_prefix="dengue"
    )

    return parser

dropped_strains = [
    'DENV1/VIETNAM/BIDV992/2006', 'DENV1/FRANCE/00475/2008', 'DENV1/VIETNAM/BIDV3990/2008', 'DENV2/HAITI/DENGUEVIRUS2HOMOSAPIENS1/2016', # Probable recombinants
    'DENV2/AUSTRALIA/QML22/2015' # Suspiciously far diverged
]

references = {
    "denv1": {"metadata": {'strain': "DENV1/NAURUISLAND/REFERENCE/1997", "accession": "NC_001477", "date": "1997-XX-XX", 'host': "NA", 'country': "Nauru", 'region': "oceania"}},
    "denv2": {"metadata": {'strain': "DENV2/THAILAND/REFERENCE/1964", "accession": "NC_001474", "date": "1964-XX-XX", 'host': "NA", 'country': "Thailand", "region": "southeast_asia"}},
    "denv3": {"metadata": {'strain': "DENV3/SRI_LANKA/REFERENCE/2000", "accession": "NC_001475", "date": "2000-XX-XX", 'host': "NA", 'country': "Sri Lanka", "region": "south_asia"}},
    "denv4": {"metadata": {'strain': "DENV4/NA/REFERENCE/2003", "accession": "NC_002640", "date": "2003-XX-XX", 'host': "NA", 'country': "NA"}},
}

for key in references:
    references[key]["genes"] = ['C', 'M', 'E', 'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', '2K', 'NS4B', 'NS5']
    references[key]["path"] = "metadata/dengue_{}_outgroup.gb".format(key)
    references[key]["include"] = 1 # include ref seq in subsampling pool
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

    config = {
        "dir": "dengue",
        "lineage": serotype,
        "file_prefix": file_prefix,
        "input_paths": None,

        "header_fields": {0:'strain', 1:'accession', 2:'date', 3:'region', 4:'country',
                        5:'division', 6: 'location', 7: 'authors', 8: 'url'},
        "filters": ( ("Dropped Strains", lambda s: s.id not in [fix_names(x) for x in dropped_strains]),
            # ("Sequence Length", lambda s: len(s.seq)>=5000),
            ("Bad Region", lambda s: s.attributes['region'] not in ['', '?']) ),

        ### Make subsampling serotype specific?? Probably not?
        "subsample": {
            "category": lambda x:(x.attributes['region'], x.attributes['date'].year, x.attributes['date'].month),
            "threshold": 3 },

        "add_urls": {
            "prefix": "https://www.ncbi.nlm.nih.gov/nuccore/%s",
            "attr": "accession" },

        "colors": ["region"],
        "color_defs": ["./colors.tsv"],
        "lat_longs": ["region"],
        "lat_long_defs": '../../fauna/source-data/geo_lat_long.tsv',
        "reference": references[serotype],
    }

    if params.sequences is not None:
        config["input_paths"] = params.sequences
    elif os.path.isfile("../../fauna/data/dengue_%s.fasta"%serotype): # Look for a serotype-specific fasta
        config["input_paths"] = ["../../fauna/data/dengue_%s.fasta"%serotype]
    else: # If it doesn't exist, try to pull serotype-specific sequences out of the all-serotype fasta (warn the user of this behavior)
        config["input_paths"] = select_serotype('../fauna/data/dengue_all.fasta', '../fauna/data/', serotype)
        print('WARNING: Did not find serotype-specific fasta file.\nPulled sequences with serotype %s from all-serotype fasta file %s\nWrote these to file %s'%(serotype, '../fauna/data/dengue.fasta', config["input_paths"]))

    years_back = params.years_back
    config['time_interval'] = [datetime.today().date(), (datetime.today()  - timedelta(days=365.25 * years_back)).date()]

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
        runner.ensure_all_segments()
        runner.subsample()
        runner.colors()
        runner.latlongs()
        runner.write_to_json()
