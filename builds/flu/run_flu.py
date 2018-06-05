import os

def run_live(
    lineages = None, resolutions = None,
    system="local",
    frequencies="complete",
    process_segment="ha",
    no_prepare=False
    ):
    lineages = ['h3n2', 'h1n1pdm', 'vic', 'yam'] if lineages is None else lineages
    resolutions = ['2y', '3y', '6y', '12y'] if resolutions is None else resolutions
    segments = ['ha', 'na', 'mp', 'np', 'ns', 'pa', 'pb1', 'pb2']

    assert process_segment in segments, "ERROR: %s is not a recognized segment."%(process_segment)

    for lineage in lineages:
        seq_files = " ".join(['../../../fauna/data/%s_%s.fasta'%(lineage, segment)
                              for segment in segments])
        for resolution in resolutions:

            if process_segment=="ha" and not no_prepare:
                call = ['python',
                    'flu.prepare.py',
                    '--lineage', lineage,
                    '--resolution', resolution,
                    '--segments', " ".join(segments),
                    '--sequences', seq_files,
                    '--titers', '../../../fauna/data/%s_public_hi_cell_titers.tsv'%(lineage),
                    '--file_prefix', 'flu_seasonal_%s_*segment*_%s'%(lineage, resolution)]
                if frequencies == "complete":
                    call = call + ['--complete_frequencies']
                print(' '.join(call))
                os.system(' '.join(call))

            call = [
                'flu.process.py',
                '--json', 'prepared/flu_seasonal_%s_ha_%s.json'%(lineage, resolution)
            ]
            if lineage == "h3n2":
                call = call + ['--annotate_fitness', '--predictors', 'cTiterSub', '--predictors_params', '1.13', '--predictors_sds', '0.72']

            if process_segment != "ha":
                import pdb; pdb.set_trace()
                call = [
                    'flu.process.py',
                    '--json', 'prepared/flu_seasonal_%s_%s_%s.json'%(lineage, process_segment, resolution)
                ]

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

def run_who(
    builds = None, lineages = None, resolutions = None,
    system="local",
    frequencies="complete",
    process_na=False,
    no_prepare = False
    ):
    builds = ['cdc', 'crick', 'niid', 'vidrl', 'who'] if builds is None else builds
    lineages = ['h3n2', 'h1n1pdm', 'vic', 'yam'] if lineages is None else lineages
    resolutions = ['2y', '6y'] if resolutions is None else resolutions
    segments = ['ha', 'na']

    for build in builds:
        for lineage in lineages:
            seq_files = " ".join(['../../../fauna/data/%s_%s.fasta'%(lineage, segment)
                                  for segment in segments])
            for resolution in resolutions:
                for passage in ['cell', 'egg']:
                    for assay in ['hi', 'fra']:

                        if lineage!='h3n2' and assay=='fra':
                            continue

                        if not (process_na or no_prepare):
                            call = ['python',
                                'flu.prepare.py',
                                '--lineage', lineage,
                                '--resolution', resolution,
                                '--segments', " ".join(segments),
                                '--sequences', seq_files,
                                '--titers', '../../../fauna/data/%s_%s_%s_%s_titers.tsv'%(lineage, build, assay, passage),
                                '--file_prefix', 'flu_%s_%s_*segment*_%s_%s_%s'%(build, lineage, resolution, passage, assay)]
                            if frequencies == "complete":
                                call = call + ['--complete_frequencies']
                            print(' '.join(call))
                            os.system(' '.join(call))

                        call = [
                            'flu.process.py',
                            '--json', 'prepared/flu_%s_%s_ha_%s_%s_%s.json'%(build, lineage, resolution, passage, assay),
                            '--titers_export'
                        ]

                        if process_na:
                            call = [
                                'flu.process.py',
                                '--json', 'prepared/flu_%s_%s_na_%s_%s_%s.json'%(build, lineage, resolution, passage, assay),
                                '--titers_export'
                            ]

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
    parser.add_argument('-v', '--version', type = str, default = 'live', help='version to run, live or who')
    parser.add_argument('-s', '--system', type = str, default = 'local', help='where to run, local, qsub or sbatch')
    parser.add_argument('-b', '--builds', nargs='+', type = str,  help ="flu builds to include")
    parser.add_argument('-l', '--lineages', nargs='+', type = str,  help ="flu lineages to include")
    parser.add_argument('-r', '--resolutions', nargs='+', type = str,  help ="flu resolutions to include")
    parser.add_argument('--frequencies', type = str, default = 'complete', help='frequencies to complete, complete or subsampled')
    parser.add_argument('--segment', type=str, default='ha',  help = "choose a segment to process; default is ha")
    parser.add_argument('--no_prepare', action="store_true", default=False,  help = "rerun previously prepared jsons")
    params = parser.parse_args()

    # only applicable to who version
    if params.builds is None:
        params.builds = ['cdc', 'crick', 'niid', 'vidrl', 'who']

    if params.lineages is None:
        params.lineages = ['h3n2', 'h1n1pdm', 'vic', 'yam']

    if params.resolutions is None:
        if params.version == "live":
            params.resolutions = ['2y', '3y', '6y', '12y']
        elif params.version == "who":
            params.resolutions = ['2y', '6y']

    if params.version == "live":
        run_live(
            lineages = params.lineages,
            resolutions = params.resolutions,
            system = params.system,
            frequencies = params.frequencies,
            process_segment = params.segment,
            no_prepare = params.no_prepare)
    elif params.version == "who":
        if params.segment == "na":
            na = True
        else:
            na = False
        run_who(
            builds = params.builds,
            lineages = params.lineages,
            resolutions = params.resolutions,
            system = params.system,
            frequencies = params.frequencies,
            process_na = na,
            no_prepare = params.no_prepare)
