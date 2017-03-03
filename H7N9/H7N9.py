from __future__ import division, print_function
from collections import defaultdict
import sys
sys.path.append('')  # need to import from base
# sys.path.append('/home/richard/Projects')  # need to import from base
from base.io_util import make_dir, remove_dir, tree_to_json, write_json, myopen
from base.sequences import sequence_set, num_date
from base.tree import tree
from base.process import process
from base.frequencies import alignment_frequencies, tree_frequencies
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import numpy as np
from datetime import datetime
from pdb import set_trace
from pprint import pprint

"""
This script is modelled on zika.py, but attempts to be more modular
"""

# HARCODED (NOT SO GOOD)
country_cmap = [
    ['?',                  "#4580CA"],
    ['canada',             "#4272CE"],
    ['china',              "#4B26B1"],
    ['czech_republic',     "#3F52CD"],
    ['guatemala',          "#5DA8A3"],
    ['hong_kong',          "#79B77D"],
    ['japan',              "#BDBB48"],
    ['spain',              "#D9AD3D"],
    ['sweden',             "#E68B35"],
    ['taiwan',             "#E56A2F"],
    ['usa',                "#DC2F24"]
]

if __name__=="__main__":

    # PARAMETERS:
    params = {
	"lineage": "H7N9",
	"input_data": "../fauna/data/H7N9",
	"proteins": ['HA'],
	# >A/chicken/Jiangxi/10877/2014|    h7n9|   EPI594171|      2014-02-18|         chicken|    china|      china|      china|      china|          egg|            other_database_import
	#        0                           1         2               3                   4           5           6           7            8          9                   10
	#       'strain',                  'virus', 'accession', 'collection_date',    'host',     'region',   'country', 'division', 'location', 'passage_category', 'submitting_lab']
	"fasta_fields": {
	    0:'strain', 2:'accession', 3:'date', 4:'host', 6:'country'
	},
	"dropped_strains": [],
	"viruses_per_month": 60,
	"geo_inference": ['country'],
	"confidence": False,
	"min_seq_length": 1500,
	"earliest_sample": datetime(2013,1,1).date()
    }

    # AUSPICE STUFF
    auspice = {
	"panels": ['tree', 'map', 'entropy'],
	"controls": {'geographic location':['country'], 'authors':['authors']},
	"color_options": {
	    "country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete", "color_map": country_cmap},
	    "num_date":{"key":"num_date", "legendTitle":"Sampling date", "menuItem":"date", "type":"continuous"},
	    "gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"}
	}
    }


    flu = process(input_data_path = params["input_data"],
		   store_data_path = 'store/' + params["lineage"] + '_',
		   build_data_path = 'build/' + params["lineage"] + '_',
		   reference='H7N9/reference_segments/4.HA.gb',
		   lat_long_fname='../fauna/source-data/geo_lat_long.tsv',
		   proteins=params['proteins'],
		   method='SLSQP', verbose=0)

    flu.load_sequences(fields=params["fasta_fields"])

    # let's find out what countries we have...
    pprint(set([x.attributes["country"] for x in flu.seqs.all_seqs.values()]))

    flu.seqs.filter(lambda s: s.attributes["country"] != "?")

    flu.seqs.filter(lambda s: s.attributes['date']>=params["earliest_sample"])
    flu.seqs.filter(lambda s: len(s.seq)>=params["min_seq_length"])
    flu.seqs.filter(lambda s: s.id not in params["dropped_strains"])
    flu.seqs.subsample(category = lambda x:(x.attributes['date'].year, x.attributes['date'].month),
	threshold=params["viruses_per_month"])
    flu.align()
    flu.build_tree()
    flu.clock_filter(n_iqd=3, plot=True)
    flu.annotate_tree(Tc=0.02, timetree=True, reroot='best', confidence=params["confidence"])
    for geo_attr in params["geo_inference"]:
	print("running geo inference for {}".format(geo_attr))
	flu.tree.geo_inference(geo_attr)
    flu.export(controls = auspice["controls"],
	       geo_attributes = params["geo_inference"],
	       color_options=auspice["color_options"],
	       panels=auspice["panels"])
