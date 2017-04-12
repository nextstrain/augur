#!/usr/bin/evn python
import os, sys, argparse, subprocess

def live_flu_build(lineage, resolution, path):
    '''
    run build for a single Virus x Resolution combination for deployment to live nextflu
    '''

    print '\n--------------------\n'
    print 'Processing lineage',lineage,'with resolution',resolution,'for all HI data'
    process = 'flu/seasonal_flu.py'
    call = map(str, [params.bin, process, '--lineage', lineage, '--resolution', resolution, '--HI', 'all_hi'])
    print call
    subprocess.call(call)

    json_types = ['entropy', 'frequencies', 'meta', 'sequences', 'tree']
    for jtype in json_types:
        f = "build/%s_all_hi_%s_%s.json"%(lineage, resolution, jtype)
        t = "%s/%s_%s_%s.json"%(path, lineage, resolution, jtype)
        os.rename(f,t)


def cdc_build(lineage, resolution, assay, passage, path):
    '''
    run build for a single Virus x Resolution combination for deployment to cdc nextflu: HI
    '''

    print '\n--------------------\n'
    print 'Processing lineage',lineage,'with resolution',resolution,'for CDC %s assays passaged in %s.'%(assay, passage)
    process = 'flu/seasonal_flu.py'
    call = map(str, [params.bin, process, '--lineage', lineage, '--resolution', resolution, '--HI cdc_%s_%s'%(assay, passage)])
    print call
    subprocess.call(call)

    json_types = ['entropy', 'frequencies', 'meta', 'sequences', 'tree', 'titer_subs_model', 'titer_tree_model', 'titers']
    for jtype in json_types:
        f = "build/%s_cdc_hi_%s_%s.json"%(lineage, resolution, jtype)
        t = "%s/%s_%s_%s.json"%(path, lineage, resolution, jtype)
        os.rename(f,t)

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

    json_types = ['entropy', 'frequencies', 'meta', 'sequences', 'tree', 'titer_subs_model', 'titer_tree_model', 'titers']
    for jtype in json_types:
        f = "build/%s_cdc_fra_%s_%s.json"%(lineage, resolution, jtype)
        t = "%s/%s_%s_%s.json"%(path, lineage, resolution, jtype)
        os.rename(f,t)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "download and process")
    parser.add_argument('--bin', type = str, default = "python")
    parser.add_argument('--builds', nargs='+', type = str, help = "builds to include")
    parser.add_argument('--lineages', nargs='+', type = str,  help = "lineages to include")
    parser.add_argument('--resolutions', nargs='+', type = str,  help = "resolutions to include")
    parser.add_argument('--assays', nargs='+', type = str, help = "assays to include")
    parser.add_argument('--passages', nargs='+', type = str, help = "passages to include")
    parser.add_argument('--live_auspice_path', default = '../../blab/nextflu/auspice/data')
    parser.add_argument('--cdc_auspice_path', default = '../../blab/nextflu-cdc/auspice/data')
    params = parser.parse_args()

    if params.builds is None:
        params.builds = ['live', 'cdc']

    if params.lineages is None:
        params.lineages = ['h3n2', 'h1n1pdm', 'vic', 'yam']

    if params.resolutions is None:
        params.resolutions = ['2y', '3y', '6y', '12y']

    if params.assays is None:
        params.assays = ['hi', 'fra']

    if params.passages is None:
        params.passages = ['egg', 'cell']

    if 'live' in params.builds:
        for lineage in params.lineages:
            for resolution in params.resolutions:
                live_flu_build(lineage, resolution, params.live_auspice_path)

    if 'cdc' in params.builds:
        for lineage in params.lineages:
            for resolutiion in params.resolutions:
                for assay in assays:
                    for passage in passages:
                        cdc_build(lineage, resolution, assay, passage, params.cdc_auspice_path)
