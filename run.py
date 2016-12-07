#!/usr/bin/env python
import os, re, time, datetime, csv, sys, argparse, subprocess

def flu_build(lineage, resolution):
	'''
	run build for single virus / resolution combination
	'''
	print '\n------------------------------\n'
	print 'Processing lineage',lineage,'with resolution',resolution
	process = 'flu/seasonal_flu.py'
	call = map(str, [params.bin, process, '--lineage', lineage, '--resolution', resolution])
	print call
	subprocess.call(call)

if __name__=="__main__":
	parser = argparse.ArgumentParser(description = "download and process")
	parser.add_argument('--bin', type = str, default = "python")
	parser.add_argument('--virus', type = str, default = "flu")
	# parser.add_argument('--zika_lineages', nargs='+', type = str,  help ="zika lineages to include")
	# parser.add_argument('--zika_resolutions', nargs='+', type = str,  help ="zika resolutions to include")
	parser.add_argument('--flu_lineages', nargs='+', type = str,  help ="flu lineages to include")
	parser.add_argument('--flu_resolutions', nargs='+', type = str,  help ="flu resolutions to include")
	params = parser.parse_args()

	# if params.zika_lineages is None:
	# 	params.zika_lineages = ['Zika']
	#
	# if params.zika_resolutions is None:
	# 	params.zika_resolutions = [None]

	if params.flu_lineages is None:
		params.flu_lineages = ['H3N2', 'H1N1pdm', 'Vic', 'Yam']

	if params.flu_resolutions is None:
		params.flu_resolutions = ['2y', '3y', '6y', '12y']

	if params.virus is "flu":
		for lineage in params.flu_lineages:
			for resolution in params.flu_resolutions:
				flu_build(lineage, resolution)
