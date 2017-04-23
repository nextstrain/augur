#!/usr/bin/evn python
import os, sys, argparse, subprocess, shutil

def live_build(lineage, resolution, path):
    '''
    run build for a single Virus x Resolution combination for deployment to live nextflu
    '''

    print '\n--------------------\n'
    print 'Processing lineage',lineage,'with resolution',resolution,'for all HI data'
    process = 'flu/seasonal_flu.py'
    call = map(str, [params.bin, process, '--lineage', lineage, '--resolution', resolution])
    print " ".join(call)
    subprocess.call(call)

    json_types = ['entropy', 'frequencies', 'sequences', 'tree', 'meta']
    for jtype in json_types:
        f = "build/%s_%s_%s.json"%(lineage, resolution, jtype)
        t = "%s%s_%s_%s.json"%(path, lineage, resolution, jtype)
        print "Copying from: ", f
        print "to: ", t
        shutil.copyfile(f,t)


def cdc_build(lineage, resolution, assay, passage, path):
    '''
    run build for a single Virus x Resolution combination for deployment to cdc nextflu: HI
    '''

    print '\n--------------------\n'
    print 'Processing lineage',lineage,'with resolution',resolution,'for CDC %s assays passaged in %s.'%(assay, passage)
    process = 'flu/seasonal_flu.py'
    call = map(str, [params.bin, process, '--lineage', lineage, '--resolution', resolution, '--HI', 'cdc_%s_%s'%(assay, passage)])
    print call
    subprocess.call(call)

    json_types = ['entropy', 'frequencies', 'meta', 'sequences', 'tree', 'titer_subs_model', 'titer_tree_model', 'titers']
    for jtype in json_types:
        f = "build/%s_cdc_%s_%s_%s_%s.json"%(lineage, assay, passage, resolution, jtype)
        t = "%s%s_cdc_%s_%s_%s_%s.json"%(path, lineage, assay, passage, resolution, jtype)
        print "Copying from: ", f
        print "to: ", t
        shutil.copyfile(f,t)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "download and process")
    parser.add_argument('--bin', type = str, default = "python")
    parser.add_argument('--builds', default=['live', 'cdc'], nargs='+', type = str, help = "builds to include")
    parser.add_argument('--lineages', default=['h3n2', 'h1n1pdm', 'vic', 'yam'], nargs='+', type = str,  help = "lineages to include")
    parser.add_argument('--resolutions', default=['2y', '3y', '6y', '12y'], nargs='+', type = str,  help = "resolutions to include")
    parser.add_argument('--assays', default=['hi', 'fra'], nargs='+', type = str, help = "assays to include")
    parser.add_argument('--passages', default=['egg', 'cell'], nargs='+', type = str, help = "passages to include")
    parser.add_argument('--live_auspice_path', default = '../nextflu/auspice/data/')
    parser.add_argument('--cdc_auspice_path', default = '../nextflu-cdc/auspice/data/')
    params = parser.parse_args()

    if 'live' in params.builds:
        for lineage in params.lineages:
            for resolution in params.resolutions:
                live_build(lineage, resolution, params.live_auspice_path)

    if 'cdc' in params.builds:
        for lineage in params.lineages:
            for resolution in params.resolutions:
                live_build(lineage, resolution, params.cdc_auspice_path)
                for assay in params.assays:
                    for passage in params.passages:
                        cdc_build(lineage, resolution, assay, passage, params.cdc_auspice_path)
