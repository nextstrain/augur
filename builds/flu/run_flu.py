import os

def build_live(
    lineages = None, resolutions = None, segments=None,
    system="local",
    frequencies="complete",
    ):
    lineages = ['h3n2', 'h1n1pdm', 'vic', 'yam'] if lineages is None else lineages
    resolutions = ['2y', '3y', '6y'] if resolutions is None else resolutions
    segments = ["ha"] if segments is None else segments

    for lineage in lineages:
        seq_files = " ".join(['../../../fauna/data/%s_%s.fasta'%(lineage, segment)
                              for segment in segments])
        for resolution in resolutions:
            call = ['python',
                'flu.prepare.py',
                '--lineage', lineage,
                '--segment', " ".join(segments),
                '--resolution', resolution,
                #'--ensure_all_segments',
                '--sequences', seq_files,
                '--titers', '../../../fauna/data/%s_hi_titers.tsv'%(lineage),
                '--file_prefix', 'flu_%s_%s'%(lineage, resolution)]
            if frequencies == "complete":
                call = call + ['--complete_frequencies']
            print(' '.join(call))
            os.system(' '.join(call))

            for segment in segments:
                call = [
                    'flu.process.py',
                    '--json', 'prepared/flu_%s_%s.json'%(lineage, resolution)]
                if (system == "qsub"):
                    call = ['qsub', 'submit_script.sh'] + call
                elif (system == "sbatch"):
                    call = ['sbatch', 'submit_flu.sh'] + call
                    # call = ['python'] + call
                    # concat = '"' + ' '.join(call) + '"'
                    # call = ['sbatch', '-n', '1', '-c', '2', '--mem', '16192', '--time', '05:59:00', '--wrap', concat]
                elif (system == "local"):
                    call = ['python'] + call
                print(' '.join(call))
                os.system(' '.join(call))

def build_cdc(
    lineages = None, resolutions = None, segments=None,
    system="local",
    frequencies="complete",
    ):
    lineages = ['h3n2', 'h1n1pdm', 'vic', 'yam'] if lineages is None else lineages
    resolutions = ['2y', '3y', '6y'] if resolutions is None else resolutions
    segments = ["ha"] if segments is None else segments

    for lineage in lineages:
        seq_files = " ".join(['../../../fauna/data/%s_%s.fasta'%(lineage, segment)
                              for segment in segments])
        for resolution in resolutions:
            for passage in ['cell', 'egg']:
                for assay in ['hi', 'fra']:
                    if lineage!='h3n2' and assay=='fra':
                        continue

                    call = ['python',
                        'flu.prepare.py',
                        '--lineage', lineage,
                        '--segment', " ".join(segments),
                        '--resolution', resolution,
                        #'--ensure_all_segments',
                        '--sequences', seq_files,
                        '--titers', '../../../fauna/data/%s_cdc_%s_%s_titers.tsv'%(lineage, assay, passage),
                        '--file_prefix', 'flu_%s_%s_%s_%s'%(lineage, resolution, passage, assay)]
                    if frequencies == "complete": # (and passage=='cell' and (titer=='hi' or lineage!='h3n2')):
                        call = call + ['--complete_frequencies']
                    print(' '.join(call))
                    os.system(' '.join(call))

                    for segment in segments:
                        call = [
                            'flu.process.py',
                            '--json', 'prepared/flu_%s_%s_%s_%s.json'%(lineage, resolution, passage, assay),
                            '--titers_export']

                        if (system == "qsub"):
                            call = ['qsub', 'submit_script.sh'] + call
                        elif (system == "sbatch"):
                            #call = ['python'] + call
                            #concat = '"' + ' '.join(call) + '"'
                            #call = ['sbatch', '-n', '1', '-c', '2', '--mem', '16192', '--time', '05:59:00', '--wrap', concat]
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
    parser.add_argument('--segments', nargs='+', type = str,  help ="flu segment to include")
    parser.add_argument('--frequencies', type = str, default = 'complete', help='frequencies to complete, complete or subsampled')
    params = parser.parse_args()

    if params.lineages is None:
        params.lineages = ['h3n2', 'h1n1pdm', 'vic', 'yam']

    if params.resolutions is None:
        if params.build == "live":
            params.resolutions = ['2y', '3y', '6y', '12y']
        elif params.build == "cdc":
            params.resolutions = ['2y', '3y', '6y']

    if params.segments is None:
        params.segments = ['ha']

    if params.build == "live":
        build_live(
            lineages = params.lineages,
            resolutions = params.resolutions,
            segments = params.segments,
            system = params.system,
            frequencies = params.frequencies)
    elif params.build == "cdc":
        build_cdc(
            lineages = params.lineages,
            resolutions = params.resolutions,
            segments = params.segments,
            system = params.system,
            frequencies = params.frequencies)
