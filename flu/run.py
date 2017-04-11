#!/usr/bin/evn python
import os, sys, argparse, subprocess

def live_flu_build(lineage, resolution, path):
    '''
    run build for a single Virus x Resolution combination for deployment to live nextflu
    '''

    print '\n--------------------\n'
    print 'Processing lineage',lineage,'with resolution',resolution,'for all HI data'
    process = 'flu/seasonal_flu.py'
    call = map(str, [params.bin, process, '--lineage', lineage, '--resolution', resolution, '--HI all_hi'])
    print call
    subprocess.call(call)

    call = "cp build/%s_all_hi_%s_entropy.json %s/%s_%s_entropy.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_all_hi_%s_frequencies.json %s/%s_%s_frequencies.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_all_hi_%s_meta.json %s/%s_%s_meta.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_all_hi_%s_sequences.json %s/%s_%s_sequences.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_all_hi_%s_tree.json %s/%s_%s_tree.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])


def cdc_hi_build(lineage, resolution, path):
    '''
    run build for a single Virus x Resolution combination for deployment to cdc nextflu: HI
    '''

    print '\n--------------------\n'
    print 'Processing lineage',lineage,'with resolution',resolution,'for CDC HI assays.'
    process = 'flu/seasonal_flu.py'
    call = map(str, [params.bin, process, '--lineage', lineage, '--resolution', resolution, '--HI cdc_hi'])
    print call
    subprocess.call(call)

    call = "cp build/%s_cdc_hi_%s_entropy.json %s/%s_%s_entropy.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_cdc_hi_%s_frequencies.json %s/%s_%s_frequencies.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_cdc_hi_%s_meta.json %s/%s_%s_meta.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_cdc_hi_%s_sequences.json %s/%s_%s_sequences.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_cdc_hi_%s_tree.json %s/%s_%s_tree.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_cdc_hi_%s_titer_subs_model.json %s/%s_%s_titer_subs_model.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_cdc_hi_%s_titer_tree_model.json %s/%s_%s_titer_tree_model.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_cdc_hi_%s_titers.json %s/%s_%s_titers.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])

def cdc_fra_build(lineage, resolution, path):
    '''
    run build for a single Virus x Resolution combination for deployment to cdc nextflu: FRA
    '''

    print '\n--------------------\n'
    print 'Processing lineage',lineage,'with resolution',resolution,'for CDC FRA assays.'
    process = 'flu/seasonal_flu.py'
    call = map(str, [params.bin, process, '--lineage', lineage, '--resolution', resolution, '--HI cdc_fra'])
    print call
    subprocess.call(call)

    call = "cp build/%s_cdc_fra_%s_entropy.json %s/%s_%s_entropy.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_cdc_fra_%s_frequencies.json %s/%s_%s_frequencies.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_cdc_fra_%s_meta.json %s/%s_%s_meta.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_cdc_fra_%s_sequences.json %s/%s_%s_sequences.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_cdc_fra_%s_tree.json %s/%s_%s_tree.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_cdc_fra_%s_titer_subs_model.json %s/%s_%s_titer_subs_model.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_cdc_fra_%s_titer_tree_model.json %s/%s_%s_titer_tree_model.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])
    call = "cp build/%s_cdc_fra_%s_titers.json %s/%s_%s_titers.json"%(lineage, resolution, path, lineage, resolution)
    print call
    subprocess.call([call])

if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "download and process")
	parser.add_argument('--bin', type = str, default = "python")
	parser.add_argument('--builds', nargs='+', type = str, help = "builds to include")
	parser.add_argument('--lineages', nargs='+', type = str,  help = "lineages to include")
	parser.add_argument('--resolutions', nargs='+', type = str,  help = "resolutions to include")
    parser.add_argument('--assays', nargs='+', type = str, help = "assays to include")
    parser.add_argument('--live_auspice_path', default = '../../blab/nextflu/auspice/data')
    parser.add_argument('--cdc_auspice_path', default = '../../blab/nextflu-cdc/auspice/data')
	params = parser.parse_args()

    if params.builds is None:
        params.builds = ['live', 'cdc']

    if params.lineages is None:
        params.lineages = ['H3N2', 'H1N1pdm', 'Vic', 'Yam']

    if params.resolutions is None:
        params.resolutions = ['2y', '3y', '6y', '12y']

    if params.assays is None:
        params.assays = ['hi', 'fra']

    if 'live' in params.builds:
        for lineage in params.lineages:
            for resolution in params.resolutions:
                live_flu_build(lineage, resolution, params.live_auspice_path)

    # TODO: add egg vs cell builds
    # if 'cdc' in params.builds:
    #     for lineage in params.lineages:
    #         for resolutiion in params.resolutiions:
    #             cdc_hi_build(lineage, resolution, params.cdc_auspice_path)
    #             cdc_fra_build(lineage, resolution, params.cdc_auspice_path)
