import os

def build_live(system="local", frequencies="complete"):
    for lineage in ['h3n2', 'h1n1pdm', 'vic', 'yam']:
        for resolution in ['2y', '3y', '6y', '12y']:

            call = ['python',
                'flu.prepare.py',
                '--lineage', lineage,
                '--resolution', resolution,
                '--sequences', '../../fauna/data/%s.fasta'%lineage,
                '--titers', '../../fauna/data/%s_hi_titers.tsv'%(lineage),
                '--file_prefix', 'flu_%s_ha_%s'%(lineage, resolution)]
            if frequencies == "complete":
                call = call + ['--complete_frequencies']
            print(' '.join(call))
            os.system(' '.join(call))

            call = [
                'flu.process.py',
                '--json', 'prepared/flu_%s_ha_%s.json'%(lineage, resolution)]
            if (system == "cluster"):
                call = ['qsub', 'submit_script.sh'] + call
            elif (system == "local"):
                call = ['python'] + call
            print(' '.join(call))
            os.system(' '.join(call))

def build_cdc(system="local", frequencies="complete"):
    for lineage in ['h3n2', 'h1n1pdm', 'vic', 'yam']:
        for resolution in ['2y', '3y', '6y']:
            for passage in ['cell', 'egg']:
                for assay in ['hi', 'fra']:
                    if lineage!='h3n2' and assay=='fra':
                        continue

                    call = ['python',
                        'flu.prepare.py',
                        '--lineage', lineage,
                        '--resolution', resolution,
                        '--sequences', '../../fauna/data/%s.fasta'%lineage,
                        '--titers', '../../fauna/data/%s_cdc_%s_%s_titers.tsv'%(lineage, assay, passage),
                        '--file_prefix', 'flu_%s_ha_%s_%s_%s'%(lineage, resolution, passage, assay)]
                    if frequencies == "complete":
                        call = call + ['--complete_frequencies']
                    print(' '.join(call))
                    os.system(' '.join(call))

                    call = [
                        'flu.process.py',
                        '--json', 'prepared/flu_%s_ha_%s_%s_%s.json'%(lineage, resolution, passage, assay),
                        '--titers_export']
                    if (system == "cluster"):
                        call = ['qsub', 'submit_script.sh'] + call
                    elif (system == "local"):
                        call = ['python'] + call
                    print(' '.join(call))
                    os.system(' '.join(call))

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run flu builds')
    parser.add_argument('-b', '--build', type = str, default = 'live', help='build to run, live or cdc')
    parser.add_argument('-s', '--system', type = str, default = 'local', help='where to run, local or cluster')
    parser.add_argument('--frequencies', type = str, default = 'complete', help='frequencies to complete, complete or subsampled')
    params = parser.parse_args()

    if params.build == "live":
        build_live(system = params.system, frequencies = params.frequencies)
    elif params.build == "cdc":
        build_cdc(system = params.system, frequencies = params.frequencies)
