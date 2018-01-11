import os

def build_live(
    lineages = None, resolutions = None,
    system="local",
    frequencies="complete",
    process_na=False
    ):
    lineages = ['h3n2', 'h1n1pdm', 'vic', 'yam'] if lineages is None else lineages
    resolutions = ['2y', '3y', '6y'] if resolutions is None else resolutions
    segments = ['ha', 'na']

    for lineage in lineages:
        seq_files = " ".join(['../../../fauna/data/%s_%s.fasta'%(lineage, segment)
                              for segment in segments])
        for resolution in resolutions:

            if not process_na:
                call = ['python',
                    'flu.prepare.py',
                    '--lineage', lineage,
                    '--resolution', resolution,
                    '--segments', " ".join(segments),
                    '--sequences', seq_files,
                    '--titers', '../../../fauna/data/%s_hi_titers.tsv'%(lineage),
                    '--file_prefix', 'flu_%s_*segment*_%s'%(lineage, resolution)]
                if frequencies == "complete":
                    call = call + ['--complete_frequencies']
                print(' '.join(call))
                os.system(' '.join(call))

            call = [
                'flu.process.py',
                '--json', 'prepared/flu_%s_ha_%s.json'%(lineage, resolution)]
            if process_na:
                call = [
                    'flu.process.py',
                    '--json', 'prepared/flu_%s_na_%s.json'%(lineage, resolution)]
            if (system == "qsub"):
                call = ['qsub', 'submit_script.sh'] + call
            elif (system == "rhino"):
                concat = '"' + ' '.join( ['python'] + call ) + '"'
                call = ['sbatch', '-n', '1', '-c', '2', '--mem', '16192', '--time', '12:00:00', '--wrap', concat]
            elif (system == "sbatch"):
                call = ['sbatch', 'submit_flu.sh'] + call
            elif (system == "local"):
                call = ['python'] + call
            print(' '.join(call))
            os.system(' '.join(call))

def build_cdc(
    lineages = None, resolutions = None,
    system="local",
    frequencies="complete",
    process_na=False
    ):
    lineages = ['h3n2', 'h1n1pdm', 'vic', 'yam'] if lineages is None else lineages
    resolutions = ['2y', '3y', '6y'] if resolutions is None else resolutions
    segments = ['ha', 'na']

    for lineage in lineages:
        seq_files = " ".join(['../../../fauna/data/%s_%s.fasta'%(lineage, segment)
                              for segment in segments])
        for resolution in resolutions:
            for passage in ['cell', 'egg']:
                for assay in ['hi', 'fra']:
                    if lineage!='h3n2' and assay=='fra':
                        continue

                    if not process_na:
                        call = ['python',
                            'flu.prepare.py',
                            '--lineage', lineage,
                            '--resolution', resolution,
                            '--segments', " ".join(segments),
                            '--sequences', seq_files,
                            '--titers', '../../../fauna/data/%s_cdc_%s_%s_titers.tsv'%(lineage, assay, passage),
                            '--file_prefix', 'flu_%s_*segment*_%s_%s_%s'%(lineage, resolution, passage, assay)]
                        if frequencies == "complete":
                            call = call + ['--complete_frequencies']
                        print(' '.join(call))
                        os.system(' '.join(call))

                    call = [
                        'flu.process.py',
                        '--json', 'prepared/flu_%s_ha_%s_%s_%s.json'%(lineage, resolution, passage, assay),
                        '--titers_export']
                    if process_na:
                        call = [
                            'flu.process.py',
                            '--json', 'prepared/flu_%s_na_%s_%s_%s.json'%(lineage, resolution, passage, assay),
                            '--titers_export']

                    if (system == "qsub"):
                        call = ['qsub', 'submit_script.sh'] + call
                    elif (system == "rhino"):
                        concat = '"' + ' '.join( ['python'] + call ) + '"'
                        call = ['sbatch', '-n', '1', '-c', '2', '--mem', '16192', '--time', '12:00:00', '--wrap', concat]
                    elif (system == "sbatch"):
                        call = ['sbatch', 'submit_flu.sh'] + call
                    elif (system == "local"):
                        call = ['python'] + call
                    print(' '.join(call))
                    os.system(' '.join(call))

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run flu builds')
    parser.add_argument('-b', '--build', type = str, default = 'live', help='build to run, live or cdc')
    parser.add_argument('-s', '--system', type = str, default = 'local', help='where to run, local, qsub or sbatch')
    parser.add_argument('-l', '--lineages', nargs='+', type = str,  help ="flu lineages to include")
    parser.add_argument('-r', '--resolutions', nargs='+', type = str,  help ="flu resolutions to include")
    parser.add_argument('--frequencies', type = str, default = 'complete', help='frequencies to complete, complete or subsampled')
    parser.add_argument('--process_na', action="store_true", default=False,  help = "supplemental run of na")
    params = parser.parse_args()

    if params.lineages is None:
        params.lineages = ['h3n2', 'h1n1pdm', 'vic', 'yam']

    if params.resolutions is None:
        if params.build == "live":
            params.resolutions = ['2y', '3y', '6y', '12y']
        elif params.build == "cdc":
            params.resolutions = ['2y', '3y', '6y']

    if params.build == "live":
        build_live(
            lineages = params.lineages,
            resolutions = params.resolutions,
            system = params.system,
            frequencies = params.frequencies,
            process_na = params.process_na)
    elif params.build == "cdc":
        build_cdc(
            lineages = params.lineages,
            resolutions = params.resolutions,
            system = params.system,
            frequencies = params.frequencies,
            process_na = params.process_na)
